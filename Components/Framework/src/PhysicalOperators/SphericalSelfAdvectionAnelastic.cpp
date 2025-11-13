/** 
 * @file SphericalSelfAdvectionAnelastic.cpp
 * @brief Source of the implementation of the spherical Coriolis term
 */


// System includes
//

// Project includes
//
#include "QuICC/PhysicalOperators/SphericalSelfAdvectionAnelastic.hpp"
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Physical {

   void SphericalSelfAdvectionAnelastic::set(Framework::Selector::PhysicalScalarField &rS,
                                             FieldComponents::Physical::Id compId, 
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Rho       = pF->evaluateLP(r, 0, 0); 
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);

      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.setSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);
               
               rS.subSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
               rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()  
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                )/ Rho(iR_) / Rho(iR_)).matrix(),  iR);
            
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.setSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);
               
               rS.subSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
               rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()  
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                )/ Rho(iR_) / Rho(iR_)).matrix(),  iR);
            
            }
         }

      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

              // Boussinesq part (not vanishing for dLogRho =0)
               rS.setSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.subSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.setSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.subSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);
               
               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.setSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.subSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.setSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.subSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
            }
         }
      }
   }

   void SphericalSelfAdvectionAnelastic::add(Framework::Selector::PhysicalScalarField &rS,
                                             FieldComponents::Physical::Id compId, 
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Rho       = pF->evaluateLP(r, 0, 0); 
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);

      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);
               
               rS.subSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
               rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()  
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                )/ Rho(iR_) / Rho(iR_)).matrix(),  iR);
            
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.addSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);
               
               rS.subSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
               rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()  
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                )/ Rho(iR_) / Rho(iR_)).matrix(),  iR);
            
            }
         }

      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

              // Boussinesq part (not vanishing for dLogRho =0)
               rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.subSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.addSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.subSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);
               
               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.subSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.addSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.subSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
            }
         }
      }
   }

   void SphericalSelfAdvectionAnelastic::sub(Framework::Selector::PhysicalScalarField &rS,
                                             FieldComponents::Physical::Id compId, 
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Rho       = pF->evaluateLP(r, 0, 0); 
      auto dLogRho   = pDF->evaluateLP(r, 0, 0); 

      if(compId == FieldComponents::Physical::R)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.subSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);
               
               rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
               rS.subSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()  
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                )/ Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               
               // Boussinesq part (not vanishing for dLogRho =0)
               rS.subSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);
               
               rS.addSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
               rS.subSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::PHI).slice(iR).array()  
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                )/ Rho(iR_) / Rho(iR_)).matrix(),  iR);
            }
         }
      } else if(compId == FieldComponents::Physical::THETA)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);               

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.subSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
            
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);               

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.subSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.addSlice((   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
            }
         }
      } else if(compId == FieldComponents::Physical::PHI)
      {
         if(c != 1.0)
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

               // Boussinesq part (not vanishing for dLogRho =0)
               rS.subSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.addSlice(c*(  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
            }
         } else
         {
            for(int iR = 0; iR < nR; ++iR)
            {
               iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
               
               // Boussinesq part (not vanishing for dLogRho =0)
               rS.subSlice((   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               rS.addSlice((   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                    * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                       / Rho(iR_) / Rho(iR_)).matrix(), iR);

               // Anelastic part (vanishing for dLogRho =0)
               rS.addSlice((  dLogRho(iR_) * ( v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                                * v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                                   ) / Rho(iR_) / Rho(iR_)).matrix(),  iR);
               
            }
         }
      }
   }

}
}
