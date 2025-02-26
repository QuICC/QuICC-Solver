/** 
 * @file SphericalSelfAdvectionAnelastic.cpp
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
#include "QuICC/PhysicalOperators/SphericalSelfAdvectionAnelastic.hpp"

// Project includes
//#include "QuICC/DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
//

namespace QuICC {

namespace Physical {

   // \todo the arguments in add, sub, and set are repeated. Should be a function
   // Issues with the return types. I got as far as
   // type:Eigen::MatrixBase<Derived>
   // But it probably should go in the .hpp file
   /*
   type:Eigen::MatrixBase<Derived> BoussinesqRcomp(const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                       const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,
                       const int iTh,
                       const int iR,
                       const int iR_, 
                       const MHDFloat c)
   {
      //return v.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array().cols();
      return c*(   v.comp(FieldComponents::Physical::R).profile(iTh,iR).array() 
                                       * w.comp(FieldComponents::Physical::PHI).profile(iTh,iR).array()
                                          / Rho(iR_) / Rho(iR_)
                                          ).matrix();

   }
   */

   void SphericalSelfAdvectionAnelastic::set(Framework::Selector::PhysicalScalarField &rS,
                                             FieldComponents::Physical::Id compId, 
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      pF->initParams(eqParams);
      auto Rho       = pF->evaluate(r, 0, 0); 
      //auto Rho       = pF->evaluate(r, 0, 0, eqParams); //overloaded version
      auto dLogRho   = pDF->evaluate(r, 0, 0);

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
      // *** to print the result ** //
      // std::cerr << "NL(R) = "<<rS.data()<<" \n";

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
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Rho       = pF->evaluate(r, 0, 0); 
      auto dLogRho   = pDF->evaluate(r, 0, 0);

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
      // *** to print the result ** //
      // std::cerr << "NL(R) = "<<rS.data()<<" \n";

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
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Rho       = pF->evaluate(r, 0, 0); 
      auto dLogRho   = pDF->evaluate(r, 0, 0); 

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
