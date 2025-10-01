/**
 * @file ISphericalHydroCfl.cpp
 * @brief Source of the CFL constraint wrapper in a spherical geometry
 */

// System includes
//

// Project includes
//
#include "QuICC/Diagnostics/ISphericalHydroCfl.hpp"

namespace QuICC {

namespace Diagnostics {

   ISphericalHydroCfl::ISphericalHydroCfl(const std::map<std::size_t,NonDimensional::SharedINumber>& params, const MHDFloat courant)
      : ISphericalCflWrapper(courant), mVelId(0)
   {
      this->mIsActive = true;
   }

   void ISphericalHydroCfl::defineVelocity(const std::size_t velId)
   {
      this->mVelId = velId;
      this->mFields.try_emplace(this->mVelId, nullptr);
   }

   void ISphericalHydroCfl::init(const std::vector<Array>& mesh)
   {
      if(this->mVelId == 0)
      {
         throw std::logic_error("Field used by CFL is not defined");
      }

      // Initialize the mesh
      this->initMesh(mesh);
   }

   Matrix ISphericalHydroCfl::initialCfl() const
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

   Matrix ISphericalHydroCfl::cfl() const
   {
      MHDFloat effVel; // Effective velocity

      const Array& r = this->mMeshSpacings.at(0);
      const Array& dr = this->mMeshSpacings.at(1);
      const Array& r_ll1 = this->mMeshSpacings.at(2);

      // Storage for minimum, global CFLs, radial CFL, horizontal CFL
      MHDFloat newCfl;
      Matrix cfl = Matrix::Constant(2, 2, std::numeric_limits<MHDFloat>::max());

      const auto& vel = this->mFields.at(this->mVelId);

      int nR = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();

      for(int i = 0; i < nR; ++i)
      {
         int iR = vel->res().cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(i);

         // Radial CFL
         effVel = vel->one().slice(i).array().abs().maxCoeff();
         newCfl = dr(iR)/effVel;
         if(newCfl < cfl(0,0))
         {
            cfl(0,0) = newCfl;
            cfl(1,0) = r(iR);
         }

         // Horizontal CFL
         effVel = (vel->two().slice(i).array().pow(2) + vel->three().slice(i).array().pow(2)).array().sqrt().maxCoeff();
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
