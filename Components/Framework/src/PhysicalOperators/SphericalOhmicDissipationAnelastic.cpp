/** 
 * @file SphericalOhmicDissipationAnelastic.cpp
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
#include "QuICC/PhysicalOperators/SphericalOhmicDissipationAnelastic.hpp"

// Project includes
//#include "QuICC/DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "DenseSM/Chebyshev/LinearMap/RadialTorPolFunction.hpp"
#include "Types/Internal/Typedefs.hpp"
//

namespace QuICC {

namespace Physical {

   void SphericalOhmicDissipationAnelastic::add(Framework::Selector::PhysicalScalarField &rS,
                                                const Resolution& res, 
                                                const Array& r, 
                                                const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                                std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for magnetic diffusivity, Eta
                                                const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                                const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Eta       = pF->evaluate(r, 0, 0); 

      for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);

            // Boussinesq part (not vanishing for dLogEta =0)
            rS.addSlice(c*Eta(iR_)*(   w.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    ).matrix(), iR);
            
            rS.addSlice(c*Eta(iR_)*(   w.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    ).matrix(), iR);

            rS.addSlice(c*Eta(iR_)*(   w.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    ).matrix(), iR);
         
         }
   }

   void SphericalOhmicDissipationAnelastic::test(Framework::Selector::PhysicalScalarField &rS,
                                                const Resolution& res, 
                                                const Array& r, 
                                                const Array& thGrid,    // Add theta grid
                                                const Array& phGrid,    // Add phi grid 
                                                const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &w,  
                                                std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for magnetic diffusivity, Eta
                                                const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                                const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto Eta       = pF->evaluate(r, 0, 0); 

      for(int iR = 0; iR < nR; ++iR)
         {
            iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
            int nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 

            // Boussinesq part (not vanishing for dLogEta =0)
            rS.setSlice(c*Eta(iR_)*(   w.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array()
                                    ).matrix(), iR);
            
            rS.addSlice(c*Eta(iR_)*(   w.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
                                    ).matrix(), iR);

            rS.addSlice(c*Eta(iR_)*(   w.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()
                                    ).matrix(), iR);
         
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
            std::cerr << "eta(iR) =  ("<<Eta(iR)<<")"<<" \n";
            std::cerr <<" \n";

            std::cerr << "j_r = "<<w.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "j_theta = "<<w.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
            std::cerr <<" \n";
            std::cerr << "j_phi = "<<w.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
            std::cerr <<" \n";

            std::cerr << "Q_j = "<< Eta(iR_)*(   w.comp(FieldComponents::Physical::R).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::R).slice(iR).array() +
                                    w.comp(FieldComponents::Physical::THETA).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::THETA).slice(iR).array() +
                              w.comp(FieldComponents::Physical::PHI).slice(iR).array() 
                                 * w.comp(FieldComponents::Physical::PHI).slice(iR).array()).matrix() <<" \n";
            std::cerr <<" \n";
            std::cerr << "c = "<< c<<" \n";
         }
   }

}
}
