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

   // Helper function to implement the calculation of Di* Q_nu/T
   // Di is the dissipation number (passed via c);
   // Q_nu = 2 nu rho (E:E -(div(v))^2/3);
   // and E is the strain rate associate to the velocity field v = u/rho
   Eigen::Matrix<MHDFloat, 
                 Eigen::Dynamic, 
                 Eigen::Dynamic>SphericalViscousDissipationAnelastic::computeViscousSlice(const int iR,
                                                                                          const int iR_,
                                                                                          const MHDFloat c,
                                                                                          const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, 
                                                                                             FieldComponents::Physical::Id>& v,
                                                                                          const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, 
                                                                                             FieldComponents::Physical::Id>& Dv,
                                                                                          const MHDFloat nu,
                                                                                          const MHDFloat T,
                                                                                          const MHDFloat Rho,
                                                                                          const MHDFloat dLogRho)
{
    return (c * 2 * Rho * nu * (
                                 // E_rr^2 
                                 ( -v.comp(FieldComponents::Physical::R).slice(iR).array()*dLogRho/Rho 
                                   + Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::R).slice(iR).array()/Rho ).pow(2)
                                 // + E_tt^2
                                 + ( Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::THETA).slice(iR).array()/Rho ).pow(2)
                                 // + E_pp^2
                                 + ( Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::PHI).slice(iR).array()/Rho ).pow(2)
                                 // + 2* (E_rt)^2
                                 + 2*( -0.5*v.comp(FieldComponents::Physical::THETA).slice(iR).array()*dLogRho/Rho  
                                      + 0.5*(  Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::THETA).slice(iR).array() 
                                                +  Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::R).slice(iR).array() )/Rho ).pow(2)
                                 // + 2* (E_rp)^2
                                 + 2*( -0.5*v.comp(FieldComponents::Physical::PHI).slice(iR).array()*dLogRho/Rho  
                                      + 0.5*(  Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::PHI).slice(iR).array() 
                                                +  Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::R).slice(iR).array() )/Rho ).pow(2)
                                 // + 2* (E_tp)^2
                                 + 2*( 0.5*(  Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::PHI).slice(iR).array() 
                                                +  Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::THETA).slice(iR).array() )/Rho ).pow(2)
                                 // -(1/3)div(v)
                                 - (1.0/3.0) * ( -v.comp(FieldComponents::Physical::R).slice(iR).array()*dLogRho/Rho ).pow(2)
                              ) / T).matrix();
}


   void SphericalViscousDissipationAnelastic::set(Framework::Selector::PhysicalScalarField &rS,
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto nu        = pV->evaluate(r, 0, 0); 
      auto T         = pT->evaluate(r, 0, 0); 
      auto Rho       = pF->evaluate(r, 0, 0); 
      auto dLogRho   = pDF->evaluate(r, 0, 0);

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);         

         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), Rho(iR_), dLogRho(iR_));

         rS.setSlice(slice, iR);


         // Test the diagonal gradient components: OK (rms is 10^-40 or so)
         /*
         rS.setSlice(((Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::R).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::THETA,FieldComponents::Physical::THETA).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::PHI,FieldComponents::Physical::PHI).slice(iR).array()
                        )).matrix(), iR);
         */
         
         // test the curl-r: NOT GOOD (rms is O(1) or so, vs O(1) rms curl values)
         // Poloidal part is ok (r curl =0)
         /*
         rS.setSlice((v.comp(FieldComponents::Physical::R).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::THETA,FieldComponents::Physical::PHI).slice(iR).array()
                        - Dv.comp(FieldComponents::Physical::PHI,FieldComponents::Physical::THETA).slice(iR).array()
                        ).matrix(), iR);
         */
         
         // test the curl-theta: OK (rms is 10^-30 or so, vs O(100) rms curl values), but check if I am cheating
         /*
         rS.setSlice(((v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        - Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::PHI).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::PHI,FieldComponents::Physical::R).slice(iR).array()
                        )).matrix(), iR);
         */
         
         // test the curl-phi: OK (rms is 10^-30 or so, vs O(100) rms curl values), but check if I am cheating
         /*
         rS.addSlice(((v.comp(FieldComponents::Physical::PHI).slice(iR).array()
                        - Dv.comp(FieldComponents::Physical::THETA,FieldComponents::Physical::R).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::THETA).slice(iR).array()
                        )).matrix(), iR);
         */
         
      }
   }

   void SphericalViscousDissipationAnelastic::add(Framework::Selector::PhysicalScalarField &rS,
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto nu        = pV->evaluate(r, 0, 0); 
      auto T         = pT->evaluate(r, 0, 0); 
      auto Rho       = pF->evaluate(r, 0, 0); 
      auto dLogRho   = pDF->evaluate(r, 0, 0);

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);         

         //rS.addSlice(c* 2*Rho(iR_)*nu(iR_) * ((
         //                                       // e_rr
         //                                       Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::R).slice(iR).array()
         //                                       )/T(iR_)).matrix(), iR);
         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), Rho(iR_), dLogRho(iR_));

         rS.addSlice(slice, iR);
      }

      
   // *** to print the result ** //
   // std::cerr << "NL(R) = "<<rS.data()<<" \n";

   
   }

   void SphericalViscousDissipationAnelastic::sub(Framework::Selector::PhysicalScalarField &rS,
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv,  
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF, // density, Rho
                                             std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pDF, // derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto nu        = pV->evaluate(r, 0, 0); 
      auto T         = pT->evaluate(r, 0, 0); 
      auto Rho       = pF->evaluate(r, 0, 0); 
      auto dLogRho   = pDF->evaluate(r, 0, 0);

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);         

         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), Rho(iR_), dLogRho(iR_));

         rS.subSlice(slice, iR);
      }
   }

}
}
