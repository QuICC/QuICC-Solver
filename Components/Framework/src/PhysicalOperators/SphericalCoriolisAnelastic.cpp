/** 
 * @file SphericalCoriolisAnelastic.cpp
 * @brief Source of the implementation of the spherical Coriolis term
 */

// Configuration includes
//

// System includes
//
// ****************
//Stuff that needs to be removed later
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <iostream>
// ****************
// External includes
//

// Class include
//
#include "QuICC/PhysicalOperators/SphericalCoriolisAnelastic.hpp"

// Project includes
//#include "QuICC/DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
//

namespace QuICC {

namespace Physical {

   void SphericalCoriolisAnelastic::set(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

                  rS.setProfile(-c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
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

                  rS.setProfile(-v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

                  rS.setProfile(-c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
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

                  rS.setProfile(-v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
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
   }

   void SphericalCoriolisAnelastic::add(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& r, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;
      int iR_;

      auto rho = pF->evaluate(r, 0, 0); 
      

      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            std::cerr << "here = " << " \n";
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(nR-1);
            std::cerr << "nR -1 = "<<nR-1 << "; iR_ ="<< iR_ << " \n";

            std::cerr << "r = "<<r(iR_) << "; vr = "<< v.comp(FieldComponents::Physical::R).profile(0,nR-1)<<" \n";

            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 

              
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

                  rS.subProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_))/rho(iR_), iTh, iR);

                  
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

                  rS.subProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

                  rS.subProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
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

                  rS.subProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
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
   }

   void SphericalCoriolisAnelastic::sub(Framework::Selector::PhysicalScalarField &rS, FieldComponents::Physical::Id compId, const Resolution& res, const Array& cosTheta, const Array& sinTheta, const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v, const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;

      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

                  rS.addProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_)), iTh, iR);
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

                  rS.addProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

                  rS.addProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_)), iTh, iR);
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

                  rS.addProfile(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*cosTheta(iTh_), iTh, iR);
               }
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
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
}
