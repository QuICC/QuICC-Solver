/**
 * @file ISphericalCflWrapper.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
 */

// Debug includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Diagnostics/ISphericalCflWrapper.hpp"

// Project includes
//
#include "QuICC/NonDimensional/CflTorsional.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/CflAlfvenDamping.hpp"
#include "QuICC/NonDimensional/CflAlfvenScale.hpp"

namespace QuICC {

namespace Diagnostics {

   ISphericalCflWrapper::ISphericalCflWrapper(const SharedIVectorWrapper spVelocity, const std::map<std::size_t,NonDimensional::SharedINumber>& params)
      : ICflWrapper(spVelocity),
        mcCourant(0.65),
        mcAlfvenScale(0),
        mcAlfvenDamping(0),
        mGlobalCfl(2)
   {
      this->mGlobalCfl.resize(params.count(NonDimensional::CflInertial::id()) + params.count(NonDimensional::CflTorsional::id()));

      int iCfl = 0;
      // Inertial wave CFL
      if(params.count(NonDimensional::CflInertial::id()) > 0)
      {
         this->mGlobalCfl(iCfl) = params.find(NonDimensional::CflInertial::id())->second->value();
         iCfl++;
      }
      // Torsional wave CFL
      if(params.count(NonDimensional::CflTorsional::id()) > 0)
      {
         this->mGlobalCfl(iCfl) = params.find(NonDimensional::CflTorsional::id())->second->value();
         iCfl++;
      }

   }

   ISphericalCflWrapper::ISphericalCflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic, const std::map<std::size_t,NonDimensional::SharedINumber>& params)
      : ICflWrapper(spVelocity, spMagnetic),
        mcCourant(0.65),
        mcAlfvenScale(params.find(NonDimensional::CflAlfvenScale::id())->second->value()),
        mcAlfvenDamping(params.find(NonDimensional::CflAlfvenDamping::id())->second->value()),
        mGlobalCfl(2)
   {
      this->mGlobalCfl.resize(params.count(NonDimensional::CflInertial::id()) + params.count(NonDimensional::CflTorsional::id()));

      int iCfl = 0;

      // Inertial wave CFL
      if(params.count(NonDimensional::CflInertial::id()) > 0)
      {
         this->mGlobalCfl(iCfl) = params.find(NonDimensional::CflInertial::id())->second->value();
         iCfl++;
      }
      // Torsional wave CFL
      if(params.count(NonDimensional::CflTorsional::id()) > 0)
      {
         this->mGlobalCfl(iCfl) = params.find(NonDimensional::CflTorsional::id())->second->value();
         iCfl++;
      }
   }

   ISphericalCflWrapper::~ISphericalCflWrapper()
   {
   }

   void ISphericalCflWrapper::init(const std::vector<Array>& mesh)
   {
      // Initialize the mesh
      this->initMesh(mesh);
   }

   void ISphericalCflWrapper::initMesh(const std::vector<Array>& mesh)
   {
      // Compute the mesh spacings
      this->mMeshSpacings.reserve(3);

      // Storage for radial grid
      this->mMeshSpacings.push_back(mesh.at(0));

      // Storage for radial grid
      this->mMeshSpacings.push_back(Array(mesh.at(0).size()));

      // Storage for horizontal average grid spacing
      this->mMeshSpacings.push_back(Array(mesh.at(0).size()));

      Array& r = this->mMeshSpacings.at(0);
      Array& dr = this->mMeshSpacings.at(1);
      Array& r_ll1 = this->mMeshSpacings.at(2);

      // Compute grid spacings
      for(int j = 0; j < r.size(); ++j)
      {
         // Get internal points grid spacing
         if(j > 0 && j < r.size() - 1)
         {
            dr(j) = std::min(std::abs(r(j) - r(j-1)), std::abs(r(j) - r(j+1)));

         // Get left endpoint grid spacing
         } else if(j > 0)
         {
            dr(j) = std::abs(r(j) - r(j-1));

            // Get right endpoint grid spacing
         } else
         {
            dr(j) = std::abs(r(j) - r(j+1));
         }

         // Compute average horizontal grid spacing
         MHDFloat effL = this->effectiveMaxL(r(j));
         r_ll1(j) = r(j)/std::sqrt(effL*(effL + 1.0));
      }
   }

   Matrix ISphericalCflWrapper::initialCfl() const
   {
      Matrix cfl = this->cfl();

      const Array& r = this->mMeshSpacings.at(0);
      const Array& dr = this->mMeshSpacings.at(1);
      const Array& r_ll1 = this->mMeshSpacings.at(2);

      // Assume a velocity of 100 to avoid problems with "zero" starting values
      MHDFloat newCfl;
      int idx;
      const int iCfl = this->mGlobalCfl.size()+1;
      newCfl = this->mcCourant*dr.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,iCfl))
      {
         cfl(0,iCfl) = newCfl;
         cfl(1,iCfl) = r(idx);
      }
      newCfl = this->mcCourant*r_ll1.minCoeff(&idx)/100.;
      if(newCfl < cfl(0,iCfl+1))
      {
         cfl(0,iCfl+1) = newCfl;
         cfl(1,iCfl+1) = r(idx);
      }

      this->updateCflMatrix(cfl);

      return cfl;
   }

   Matrix ISphericalCflWrapper::cfl() const
   {
      MHDFloat effVel; // Effective velocity

      const Array& r = this->mMeshSpacings.at(0);
      const Array& dr = this->mMeshSpacings.at(1);
      const Array& r_ll1 = this->mMeshSpacings.at(2);

      // Storage for minimum, global CFLs, radial CFL, horizontal CFL
      MHDFloat newCfl;
      Matrix cfl = Matrix::Constant(2,this->mGlobalCfl.size()+2+1, std::numeric_limits<MHDFloat>::max());
      int iCfl = this->mGlobalCfl.size()+1;

      // Set global CFL
      if(this->mGlobalCfl.size() > 0)
      {
         cfl.row(0).segment(1,this->mGlobalCfl.size()) = this->mGlobalCfl.transpose();
         cfl.row(1).segment(1,this->mGlobalCfl.size()).array() = -1.0;
      }

      int nR = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();

      if(this->mspVelocity && this->mspMagnetic)
      {
         MHDFloat aD;
         Matrix p;
         for(int i = 0; i < nR; ++i)
         {
            int iR = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(i);

            // Radial CFL
            aD = std::pow(this->mcAlfvenDamping/dr(iR),2);
            p = this->mspMagnetic->one().slice(i).array().pow(2)*this->mcAlfvenScale;
            effVel = (p.array()/(p.array() + aD).array().sqrt() + this->mspVelocity->one().slice(i).array().abs()).maxCoeff();
            newCfl = dr(iR)/effVel;
            if(newCfl < cfl(0,iCfl))
            {
               cfl(0,iCfl) = newCfl;
               cfl(1,iCfl) = r(iR);
            }

            // Horizontal CFL
            aD = std::pow(this->mcAlfvenDamping/r_ll1(iR),2);
            p = (this->mspMagnetic->two().slice(i).array().pow(2) + this->mspMagnetic->three().slice(i).array().pow(2))*this->mcAlfvenScale;
            effVel = (p.array()/(p.array() + aD).array().sqrt() + (this->mspVelocity->two().slice(i).array().pow(2) + this->mspVelocity->three().slice(i).array().pow(2)).array().sqrt()).maxCoeff();
            newCfl = r_ll1(iR)/effVel;
            if(newCfl < cfl(0,iCfl+1))
            {
               cfl(0,iCfl+1) = newCfl;
               cfl(1,iCfl+1) = r(iR);
            }
         }

      } else
      {
         for(int i = 0; i < nR; ++i)
         {
            int iR = this->mspVelocity->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(i);

            // Radial CFL
            effVel = this->mspVelocity->one().slice(i).array().abs().maxCoeff();
            newCfl = dr(iR)/effVel;
            if(newCfl < cfl(0,iCfl))
            {
               cfl(0,iCfl) = newCfl;
               cfl(1,iCfl) = r(iR);
            }

            // Horizontal CFL
            effVel = (this->mspVelocity->two().slice(i).array().pow(2) + this->mspVelocity->three().slice(i).array().pow(2)).array().sqrt().maxCoeff();
            newCfl = r_ll1(iR)/effVel;
            if(newCfl < cfl(0,iCfl+1))
            {
               cfl(0,iCfl+1) = newCfl;
               cfl(1,iCfl+1) = r(iR);
            }
         }
      }

      cfl.row(0).tail(cfl.cols()-1).array() *= this->mcCourant;
      this->updateCflMatrix(cfl);

      return cfl;
   }

}
}
