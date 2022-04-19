/** 
 * @file SphericalSComponent.cpp
 * @brief Source of the implementation of the spherical Z component of a field
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalOperators/SphericalSComponent.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

   void SphericalSComponent::set(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      if(c != 1.0)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.setProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
            }
         }
      } else
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.setProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
            }
         }
      }
   }

   void SphericalSComponent::add(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      if(c != 1.0)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.addProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               rS.addProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
            }
         }
      } else
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.addProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               rS.addProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
            }
         }
      }
   }

   void SphericalSComponent::sub(Framework::Selector::PhysicalScalarField &rS, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      if(c != 1.0)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.subProfile(c*(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
               rS.subProfile(c*(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
            }
         }
      } else
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.subProfile(v.comp(FieldComponents::Physical::R).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               rS.subProfile(v.comp(FieldComponents::Physical::THETA).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
            }
         }
      }
   }

}
}
