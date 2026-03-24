/**
 * @file ISphericalMagneticCfl.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/ISphericalMagneticCfl.hpp"
#include "QuICC/NonDimensional/CflAlfvenDamping.hpp"
#include "QuICC/NonDimensional/CflAlfvenScale.hpp"

#include "Environment/Cfl.hpp"
#ifdef QUICC_HAS_CUDA_BACKEND
#include "Cuda/CudaUtil.hpp"
#endif
namespace QuICC {

namespace Diagnostics {

   ISphericalMagneticCfl::ISphericalMagneticCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant)
      : ISphericalCflWrapper(courant),
        mcAlfvenScale((params.count(NonDimensional::CflAlfvenScale::id()) > 0) ? params.find(NonDimensional::CflAlfvenScale::id())->second->value() : 0),
        mcAlfvenDamping((params.count(NonDimensional::CflAlfvenDamping::id()) > 0) ? params.find(NonDimensional::CflAlfvenDamping::id())->second->value() : 0),
        mVelId(0),
        mMagId(0)
   {
      if(this->mcAlfvenScale == 0 || this->mcAlfvenDamping == 0)
      {
         throw std::logic_error("Alfven wave parameters are missing");
      }

      this->mIsActive = true;
   }

   void ISphericalMagneticCfl::defineVelocity(const std::size_t velId)
   {
      this->mVelId = velId;
      this->mFields.try_emplace(this->mVelId, nullptr);
   }

   void ISphericalMagneticCfl::defineMagnetic(const std::size_t magId)
   {
      this->mMagId = magId;
      this->mFields.try_emplace(this->mMagId, nullptr);
   }

   void ISphericalMagneticCfl::init(const std::vector<Array>& mesh)
   {
      if(this->mVelId == 0 || this->mMagId == 0)
      {
         throw std::logic_error("Fields used by CFL are not defined");
      }

      // Initialize the mesh
      this->initMesh(mesh);
   }

   Matrix ISphericalMagneticCfl::initialCfl() const
   {
      const Array& r = this->mMeshSpacings.at(0);
      const Array& dr = this->mMeshSpacings.at(1);
      const Array& r_ll1 = this->mMeshSpacings.at(2);

      int nR = this->mMeshSpacings.at(0).size();
      QuICC::Cfl_nR = nR;
      if (!QuICC::Cfl_r)
        cudaMalloc((void**)&QuICC::Cfl_r, nR* sizeof(double));
      if (!QuICC::Cfl_dr)
        cudaMalloc((void**)&QuICC::Cfl_dr, nR* sizeof(double));
      if (!QuICC::Cfl_r_ll1)
        cudaMalloc((void**)&QuICC::Cfl_r_ll1, nR* sizeof(double));

      if (!QuICC::newCfl_radial)
        cudaMalloc((void**)&QuICC::newCfl_radial, 2* sizeof(double));
      if (!QuICC::newCfl_horizontal)
        cudaMalloc((void**)&QuICC::newCfl_horizontal, 2* sizeof(double));
      double* temp_arr = (double*)calloc(nR, sizeof(double));
      for (int i = 0; i < nR; ++i)
      {
         temp_arr[i] = this->mMeshSpacings.at(0)(i);
      }
      cudaMemcpy(QuICC::Cfl_r, temp_arr, nR * sizeof(double),
         cudaMemcpyHostToDevice);
      for (int i = 0; i < nR; ++i)
      {
         temp_arr[i] = this->mMeshSpacings.at(1)(i);
      }
      cudaMemcpy(QuICC::Cfl_dr, temp_arr, nR * sizeof(double),
         cudaMemcpyHostToDevice);
      for (int i = 0; i < nR; ++i)
      {
         temp_arr[i] = this->mMeshSpacings.at(2)(i);
      }
      cudaMemcpy(QuICC::Cfl_r_ll1, temp_arr, nR * sizeof(double),
         cudaMemcpyHostToDevice);
      free(temp_arr);

      QuICC::Cfl_mcAlfvenDamping = this->mcAlfvenDamping;
      QuICC::Cfl_mcAlfvenScale = this->mcAlfvenScale;

      
      Matrix cfl = Matrix::Constant(2,2, std::numeric_limits<MHDFloat>::max());
      // Assume a velocity of 100 to avoid problems with "zero" starting values
      MHDFloat newCfl;
      int idx;

      newCfl = this->mcCourant*dr.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,0))
      {
         cfl(0,0) = newCfl;
         cfl(1,0) = r(idx);
      }
      newCfl = this->mcCourant*r_ll1.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,1))
      {
         cfl(0,1) = newCfl;
         cfl(1,1) = r(idx);
      }

      return cfl;
   }

   Matrix ISphericalMagneticCfl::cfl() const
   {
      MHDFloat effVel; // Effective velocity

      const Array& r = this->mMeshSpacings.at(0);
      const Array& dr = this->mMeshSpacings.at(1);
      const Array& r_ll1 = this->mMeshSpacings.at(2);

      // Storage for minimum, global CFLs, radial CFL, horizontal CFL
      MHDFloat newCfl;
      Matrix cfl = Matrix::Constant(2,2, std::numeric_limits<MHDFloat>::max());

      if (1)
      {
        double radialCfl[2];
        cudaMemcpy(radialCfl, QuICC::newCfl_radial, 2 * sizeof(double),
         cudaMemcpyDeviceToHost);
      cfl(0, 0) = radialCfl[0];
      cfl(1, 0) = radialCfl[1];

      double horizontalCfl[2];
      cudaMemcpy(horizontalCfl, QuICC::newCfl_horizontal, 2 * sizeof(double),
         cudaMemcpyDeviceToHost);

      cfl(0, 1) = horizontalCfl[0];
      cfl(1, 1) = horizontalCfl[1];

      }
   else
      {
        const auto& vel = this->mFields.at(this->mVelId);
        const auto& mag = this->mFields.at(this->mMagId);

        int nR = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();

        MHDFloat aD;
        Matrix p;
        for (int j = 0; j < 10; ++j)
        {
            //printf("cpu %d %e\n", j, vel->one().slice(0).array()(j,0));
       
       }
        for(int i = 0; i < nR; ++i)
        {
            int iR = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(i);

            // Radial CFL
            aD = std::pow(this->mcAlfvenDamping/dr(iR),2);
            p = mag->one().slice(i).array().pow(2)*this->mcAlfvenScale;
            
               
            effVel = (p.array()/(p.array() + aD).array().sqrt() + vel->one().slice(i).array().abs()).maxCoeff();
            newCfl = dr(iR)/effVel;
            if(newCfl < cfl(0,0))
            {
            cfl(0,0) = newCfl;
            cfl(1,0) = r(iR);
            }

            // Horizontal CFL
            aD = std::pow(this->mcAlfvenDamping/r_ll1(iR),2);
            p = (mag->two().slice(i).array().pow(2) + mag->three().slice(i).array().pow(2))*this->mcAlfvenScale;
            effVel = (p.array()/(p.array() + aD).array().sqrt() + (vel->two().slice(i).array().pow(2) + vel->three().slice(i).array().pow(2)).array().sqrt()).maxCoeff();
            newCfl = r_ll1(iR)/effVel;
            if(newCfl < cfl(0,1))
            {
            cfl(0,1) = newCfl;
            cfl(1,1) = r(iR);
            }
            /* for (int k = 0; k < 64; ++k)
             for (int j = 0; j < 96; ++j)
        {
                 double val = (p(j,k)/sqrt(p(j,k) + aD) + sqrt(vel->two().slice(i).array()(j,k) * vel->two().slice(i).array()(j,k)  + vel->three().slice(i).array()(j,k) * vel->three().slice(i).array()(j,k)));
            if (r_ll1(iR)/val < 6e-5)printf("cpu %d %d %e %e %e\n", j, k, r_ll1(iR), val, r_ll1(iR)/val) ;
        }*/
        }
      }
      //  printf("gpu %e %e %e %e\n", radialCfl[0], radialCfl[1], horizontalCfl[0], horizontalCfl[1]);
      //printf("corr %e %e %e %e\n", cfl(0, 0), cfl(1, 0), cfl(0, 1), cfl(1, 1));

      cfl.row(0).array() *= this->mcCourant;

      return cfl;
   }

}
}
