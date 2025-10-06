/** 
 * @file SphericalLorentzAnelastic.cpp
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
#include "QuICC/PhysicalOperators/SphericalLorentzAnelastic.hpp"

// Project includes
//#include "QuICC/DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
//

namespace QuICC {

namespace Physical {

   void SphericalLorentzAnelastic::set(Framework::Selector::PhysicalScalarField &rS,
                                             FieldComponents::Physical::Id compId, 
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Rho       = pF->evaluate(r, 0, 0); 

      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
            
            rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
         
         }


      } else if(compId == FieldComponents::Physical::THETA)
      {

         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
            
         }
         
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
         }

      }
   }


   void SphericalLorentzAnelastic::add(Framework::Selector::PhysicalScalarField &rS,
                                             FieldComponents::Physical::Id compId, 
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Rho       = pF->evaluate(r, 0, 0); 

      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.subSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
            
            rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
         
         }


      } else if(compId == FieldComponents::Physical::THETA)
      {

         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.subSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
            
         }
         
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.subSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
         }

      }
   }

   void SphericalLorentzAnelastic::test(Framework::Selector::PhysicalScalarField &rS,
                                             FieldComponents::Physical::Id compId, 
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Array& thGrid,    // Add theta grid
                                             const Array& phGrid,    // Add phi grid 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Rho       = pF->evaluate(r, 0, 0); 

      if(compId == FieldComponents::Physical::R)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            int nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
            
            rS.addSlice(c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);
            
            // to test the implementation
            std::cerr << "(iR, iR_) = ("<<iR<<","<<iR_<<")"<<" \n";
            std::cerr << "(r(iR), r(iR_)) = ("<<r(iR)<<","<<r(iR_)<<")"<<" \n";
            // Print theta and phi coordinates
            std::cerr << " theta = \n";
            
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               int iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);
               MHDFloat theta = thGrid(iTh_);
               
               std::cerr << theta << " ";
               
            }
            
            // Print phi values - use grid size directly since phi is typically uniform
            int nPh = phGrid.size();
            std::cerr << "\n phi = \n";
            for(int iPh = 0; iPh < nPh; ++iPh)
            {
               MHDFloat phi = phGrid(iPh);
               std::cerr <<phi << " ";
            }
            
            std::cerr << "\n";
            // for density_type=0, this is r
            std::cerr << "rho(iR) =  ("<<Rho(iR)<<")"<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_theta = "<<v.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_phi = "<<v.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "w_theta = "<<w.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "w_phi = "<<w.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            
            std::cerr << "rnLcomp = "<< -c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / Rho(iR_)).matrix() + c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / Rho(iR_)).matrix() <<" \n";
            
            
         }


      } else if(compId == FieldComponents::Physical::THETA)
      {

         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            int nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);

            // to test the implementation
            std::cerr << "(iR, iR_) = ("<<iR<<","<<iR_<<")"<<" \n";
            std::cerr << "(r(iR), r(iR_)) = ("<<r(iR)<<","<<r(iR_)<<")"<<" \n";
            // Print theta and phi coordinates
            std::cerr << " theta = \n";
            
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               int iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);
               MHDFloat theta = thGrid(iTh_);
               
               std::cerr << theta << " ";
               
            }
            
            // Print phi values - use grid size directly since phi is typically uniform
            int nPh = phGrid.size();
            std::cerr << "\n phi = \n";
            for(int iPh = 0; iPh < nPh; ++iPh)
            {
               MHDFloat phi = phGrid(iPh);
               std::cerr <<phi << " ";
            }
            
            std::cerr << "\n";
            // for density_type=0, this is r
            std::cerr << "rho(iR) =  ("<<Rho(iR)<<")"<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_theta = "<<v.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_phi = "<<v.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "w_theta = "<<w.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "w_phi = "<<w.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            
            std::cerr << "rnLcomp = "<< -c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    / Rho(iR_)).matrix() + c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / Rho(iR_)).matrix() <<" \n";

            
            
         }
         
      } else if(compId == FieldComponents::Physical::PHI)
      {
         for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            int nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 

            // Boussinesq part (not vanishing for dLogRho =0)
            rS.setSlice(-c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);

            rS.addSlice(c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / Rho(iR_)).matrix(), iR);

            // to test the implementation
            std::cerr << "(iR, iR_) = ("<<iR<<","<<iR_<<")"<<" \n";
            std::cerr << "(r(iR), r(iR_)) = ("<<r(iR)<<","<<r(iR_)<<")"<<" \n";
            // Print theta and phi coordinates
            std::cerr << " theta = \n";
            
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               int iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);
               MHDFloat theta = thGrid(iTh_);
               
               std::cerr << theta << " ";
               
            }
            
            // Print phi values - use grid size directly since phi is typically uniform
            int nPh = phGrid.size();
            std::cerr << "\n phi = \n";
            for(int iPh = 0; iPh < nPh; ++iPh)
            {
               MHDFloat phi = phGrid(iPh);
               std::cerr <<phi << " ";
            }
            
            std::cerr << "\n";
            // for density_type=0, this is r
            std::cerr << "rho(iR) =  ("<<Rho(iR)<<")"<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_theta = "<<v.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "v_phi = "<<v.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "w_theta = "<<w.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "w_phi = "<<w.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            
            std::cerr << "rnLcomp = "<< -c*(   v.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    / Rho(iR_)).matrix() + c*(   v.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    / Rho(iR_)).matrix() <<" \n";

            
         }

      }

   }

}
}
