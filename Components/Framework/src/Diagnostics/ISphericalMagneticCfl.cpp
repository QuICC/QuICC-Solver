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
      Matrix cfl = this->cfl();

      const Array& r = this->mMeshSpacings.at(0);
      const Array& dr = this->mMeshSpacings.at(1);
      const Array& r_ll1 = this->mMeshSpacings.at(2);

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

      const auto& vel = this->mFields.at(this->mVelId);
      const auto& mag = this->mFields.at(this->mMagId);

      int nR = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();

      MHDFloat aD;
      Matrix p;
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
      }

      cfl.row(0).array() *= this->mcCourant;

      return cfl;
   }

}
}
