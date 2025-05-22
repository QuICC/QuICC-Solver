/** 
 * @file SphericalViscousDissipationAnelastic.cpp
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
#include "QuICC/PhysicalOperators/SphericalViscousDissipationAnelastic.hpp"

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
                       const Datatypes::SymmetricTensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv,
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

   void SphericalViscousDissipationAnelastic::set(Framework::Selector::PhysicalScalarField &rS,
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Array& cosTheta, 
                                             const Array& sinTheta,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::SymmetricTensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pV,
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pT,
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   { /* TO BE IIMPLEMENTED IF NEEDED*/}

   void SphericalViscousDissipationAnelastic::add(Framework::Selector::PhysicalScalarField &rS,
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Array& cosTheta, 
                                             const Array& sinTheta,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::SymmetricTensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pV,
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pT,
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh;
      int iTh_;
      int iR_;

      auto nu        = pV->evaluate(r, 0, 0); 
      auto T         = pT->evaluate(r, 0, 0); 
      auto Rho       = pF->evaluate(r, 0, 0); 
      auto dLogRho   = pDF->evaluate(r, 0, 0);

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);
         nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 
         
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            iTh_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR);

            rS.addProfile(c* 2*Rho(iR_)*nu(iR_) * (
                                                   // e_rr = dvr/dr
                                                   Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::R).profile(iTh,iR)
                                                   )/T(iR_), iTh, iR);
            
            // coriolis example:
            // rS.subProfile(c*(v.comp(FieldComponents::Physical::PHI).profile(iTh,iR)*sinTheta(iTh_))/rho(iR_), iTh, iR);
            //advection example:
            //rS.setSlice(c*(   v.comp(FieldComponents::Physical::PHI).slice(iR).array() 
            // * w.comp(FieldComponents::Physical::THETA).slice(iR).array()
            // / Rho(iR_) / Rho(iR_)).matrix(), iR);
         }
      }
         
      
      
      
   // *** to print the result ** //
   // std::cerr << "NL(R) = "<<rS.data()<<" \n";

   
   }

   void SphericalViscousDissipationAnelastic::sub(Framework::Selector::PhysicalScalarField &rS,
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Array& cosTheta, 
                                             const Array& sinTheta,
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::SymmetricTensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pV,
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pT,
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // intended for density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // intended for derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   { /* TO BE IIMPLEMENTED IF NEEDED*/}

}
}
