/**
 * @file ISphericalMagneticInviscidCfl.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/ISphericalMagneticInviscidCfl.hpp"
#include "QuICC/NonDimensional/CflAlfvenDamping.hpp"
#include "QuICC/NonDimensional/CflAlfvenScale.hpp"

namespace QuICC {

namespace Diagnostics {

   ISphericalMagneticInviscidCfl::ISphericalMagneticInviscidCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant)
      : ISphericalInviscidCflWrapper(courant),
        mVelId(0),
        mMagId(0)
   {
      this->mIsActive = true;
   }

   void ISphericalMagneticInviscidCfl::defineVelocity(const std::size_t velId)
   {
      this->mVelId = velId;
      this->mFields.try_emplace(this->mVelId, nullptr);
   }

   void ISphericalMagneticInviscidCfl::defineMagnetic(const std::size_t magId)
   {
      this->mMagId = magId;
      this->mFields.try_emplace(this->mMagId, nullptr);
   }

   void ISphericalMagneticInviscidCfl::init(const std::vector<Array>& mesh)
   {
      if(this->mVelId == 0 || this->mMagId == 0)
      {
         throw std::logic_error("Fields used by CFL are not defined");
      }

      // Initialize the mesh
      this->initMesh(mesh);
   }

   Matrix ISphericalMagneticInviscidCfl::initialCfl() const
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

   Matrix ISphericalMagneticInviscidCfl::cfl() const
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

      for(int i = 0; i < nR; ++i)
      {
         int iR = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(i);

         // Radial CFL
         effVel = (dr(iR)/mag->one().slice(i).array().abs()).minCoeff();
         newCfl = std::pow(effVel,2);
         if(newCfl < cfl(0,0))
         {
            cfl(0,0) = newCfl;
            cfl(1,0) = r(iR);
         }

         // Horizontal CFL
         effVel = (r_ll1(iR)/(mag->two().slice(i).array().pow(2) + mag->three().slice(i).array().pow(2)).array().sqrt()).minCoeff();
         newCfl = std::pow(effVel,2);
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
