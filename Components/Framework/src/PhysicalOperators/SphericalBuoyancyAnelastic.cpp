/** 
 * @file SphericalBuoyancyAnelastic.cpp
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
#include "QuICC/PhysicalOperators/SphericalBuoyancyAnelastic.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

   void SphericalBuoyancyAnelastic::add(Framework::Selector::PhysicalScalarField &rS, 
                                        FieldComponents::Physical::Id compId, 
                                        const Resolution& res, 
                                        const Array& r, 
                                        const Framework::Selector::PhysicalScalarField &q, 
                                        std::shared_ptr<QuICC::DenseSM::IGenericProfile> pAlpha,
                                        std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT,
                                        std::shared_ptr<QuICC::DenseSM::IGenericProfile> pCp,
                                        std::shared_ptr<QuICC::DenseSM::IGenericProfile> pG,
                                        const MHDFloat c)
   {
      auto alpha   = pAlpha->evaluateLP(r, 0, 0); 
      auto temp    = pT->evaluateLP(r, 0, 0); 
      auto cp      = pCp->evaluateLP(r, 0, 0); 
      auto gravity = pG->evaluateLP(r, 0, 0); 

      if(compId == FieldComponents::Physical::R)
      {
         int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int iR_;

         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               rS.addSlice(c*alpha(iR_)*temp(iR_)*gravity(iR_)*q.slice(iR)/cp(iR_), iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               rS.addSlice(alpha(iR_)*temp(iR_)*gravity(iR_)*q.slice(iR)/cp(iR_), iR);
            }
         }
      }
   }

}
}
