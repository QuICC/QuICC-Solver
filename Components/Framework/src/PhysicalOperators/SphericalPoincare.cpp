/** 
 * @file SphericalPoincare.cpp
 * @brief Source of the implementation of the spherical Poincare term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalOperators/SphericalPoincare.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

   void SphericalPoincare::set(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         rS.setZeros();
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = -csa*rGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               rS.setProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.setProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
            }
         }
      }
   }

   void SphericalPoincare::add(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         //
         // Zero
         // 
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               rS.subProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.addProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
            }
         }
      }
   }

   void SphericalPoincare::sub(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& rGrid, const Array& thGrid, const Array& phGrid, const MHDFloat t, const MHDFloat alpha, const MHDFloat c)
   {
      MHDFloat csa = c*std::sin(alpha);
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         // 
         // Zero
         // 
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               rS.addProfile((coeff*(phGrid.array() + t).array().cos()).matrix(), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            coeff = csa*rGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

               rS.subProfile((coeff*std::cos(thGrid(iTh_))*(phGrid.array() + t).array().sin()).matrix(), iTh, iR);
            }
         }
      }
   }

}
}
