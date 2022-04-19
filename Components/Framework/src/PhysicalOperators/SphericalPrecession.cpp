/** 
 * @file SphericalPrecession.cpp
 * @brief Source of the implementation of the spherical Coriolis + precession term
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalOperators/SphericalPrecession.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

   void SphericalPrecession::set(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      MHDFloat cA = c*preC*std::sin(alpha);
      MHDFloat cB = c*(corC + preC*std::cos(alpha));
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // Theta component
               rS.setProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
               // Phi components
               coeff = cA*std::cos(theta);
               rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::sin(theta); 
               rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // R component
               rS.setProfile((-cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
               // Phi components
               coeff = cA*std::sin(theta);
               rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::cos(theta);
               rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // R components
               coeff = -cA*std::cos(theta);
               rS.setProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::sin(theta);
               rS.addProfile(coeff*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
               // Theta components
               coeff = cA*std::sin(theta);
               rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::cos(theta);
               rS.addProfile(coeff*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);
            }
         }
      }
   }

   void SphericalPrecession::add(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      MHDFloat cA = c*preC*std::sin(alpha);
      MHDFloat cB = c*(corC + preC*std::cos(alpha));
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // Theta component
               rS.addProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
               // Phi components
               coeff = cA*std::cos(theta);
               rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::sin(theta); 
               rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);

            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // R component
               rS.subProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
               // Phi components
               coeff = cA*std::sin(theta);
               rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::cos(theta);
               rS.subProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // R components
               coeff = cA*std::cos(theta);
               rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::sin(theta);
               rS.addProfile(coeff*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
               // Theta components
               coeff = cA*std::sin(theta);
               rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::cos(theta);
               rS.addProfile(coeff*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);

            }
         }
      }
   }

   void SphericalPrecession::sub(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& thGrid, const Array& phGrid, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat t, const MHDFloat alpha, const MHDFloat corC, const MHDFloat preC, const MHDFloat c)
   {
      MHDFloat cA = c*preC*std::sin(alpha);
      MHDFloat cB = c*(corC + preC*std::cos(alpha));
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;

      MHDFloat theta;
      MHDFloat coeff;
      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // Theta component
               rS.subProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
               // Phi components
               coeff = cA*std::cos(theta);
               rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::sin(theta); 
               rS.addProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // R component
               rS.addProfile((cA*(phGrid.array() + t).array().sin()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
               // Phi components
               coeff = cA*std::sin(theta);
               rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::cos(theta);
               rS.addProfile(coeff*v.comp(FieldComponents::Physical::PHI).profile(iTh,iR), iTh, iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               theta = thGrid(res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

               // R components
               coeff = cA*std::cos(theta);
               rS.addProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::R).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::sin(theta);
               rS.subProfile(coeff*v.comp(FieldComponents::Physical::R).profile(iTh,iR), iTh, iR);
               // Theta components
               coeff = cA*std::sin(theta);
               rS.subProfile((coeff*(phGrid.array() + t).array().cos()*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR).array()).matrix(), iTh, iR);
               coeff = cB*std::cos(theta);
               rS.subProfile(coeff*v.comp(FieldComponents::Physical::THETA).profile(iTh,iR), iTh, iR);
            }
         }
      }
   }

}
}
