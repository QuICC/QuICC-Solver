/** 
 * @file SphericalBuoyancy.cpp
 * @brief Source of the implementation of the spherical buoyancy term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalOperators/SphericalBuoyancy.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

   void SphericalBuoyancy::set(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const Framework::Selector::PhysicalScalarField &q, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int iR_;
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               rS.setSlice(c*r(iR_)*q.slice(iR), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               rS.setSlice(r(iR_)*q.slice(iR), iR);
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         rS.setZeros();
      } else if(compId == FieldComponents::Physical::PHI)
      {
         rS.setZeros();
      }
   }

   void SphericalBuoyancy::add(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const Framework::Selector::PhysicalScalarField &q, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int iR_;

         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               rS.addSlice(c*r(iR_)*q.slice(iR), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               rS.addSlice(r(iR_)*q.slice(iR), iR);
            }
         }
      }
   }

   void SphericalBuoyancy::sub(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const Framework::Selector::PhysicalScalarField &q, const MHDFloat c)
   {
      if(compId == FieldComponents::Physical::R)
      {
         int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int iR_;
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               rS.subSlice(c*r(iR_)*q.slice(iR), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               rS.subSlice(r(iR_)*q.slice(iR), iR);
            }
         }
      }
   }

}
}
