/** 
 * @file SphericalViscousDissipationAnelastic.cpp
 * @brief Source of the implementation of the spherical Coriolis term
 */

// System includes
//
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <iostream>

// Project includes
#include "QuICC/PhysicalOperators/SphericalViscousDissipationAnelastic.hpp"
#include "DenseSM/IGenericProfile.hpp"
#include "Types/Typedefs.hpp"
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
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // density, Rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto nu        = pV->evaluateLP(r, 0, 0); 
      auto T         = pT->evaluateLP(r, 0, 0); 
      auto Rho       = pF->evaluateLP(r, 0, 0); 
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);

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
         
         // test the curl-r: OK
         // Poloidal part is ok (r curl =0)
         /*
         rS.setSlice((v.comp(FieldComponents::Physical::R).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::THETA,FieldComponents::Physical::PHI).slice(iR).array()
                        - Dv.comp(FieldComponents::Physical::PHI,FieldComponents::Physical::THETA).slice(iR).array()
                        ).matrix(), iR);
         */
         
         // test the curl-theta: OK
         /*
         rS.setSlice(((v.comp(FieldComponents::Physical::THETA).slice(iR).array()
                        - Dv.comp(FieldComponents::Physical::R,FieldComponents::Physical::PHI).slice(iR).array()
                        + Dv.comp(FieldComponents::Physical::PHI,FieldComponents::Physical::R).slice(iR).array()
                        )).matrix(), iR);
         */         
         
         // test the curl-phi: OK?
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
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // density, Rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto nu        = pV->evaluateLP(r, 0, 0); 
      auto T         = pT->evaluateLP(r, 0, 0); 
      auto Rho       = pF->evaluateLP(r, 0, 0); 
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);         

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
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // density, Rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      auto nu        = pV->evaluateLP(r, 0, 0); 
      auto T         = pT->evaluateLP(r, 0, 0); 
      auto Rho       = pF->evaluateLP(r, 0, 0); 
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);         

         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), Rho(iR_), dLogRho(iR_));

         rS.subSlice(slice, iR);
      }
   }

   void SphericalViscousDissipationAnelastic::test(Framework::Selector::PhysicalScalarField &rS,
                                             const Resolution& res, 
                                             const Array& r, 
                                             const Array& thGrid,    // Add theta grid
                                             const Array& phGrid,    // Add phi grid 
                                             const Datatypes::VectorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &v,
                                             const Datatypes::TensorField<Framework::Selector::PhysicalScalarField, FieldComponents::Physical::Id> &Dv,  
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pV, // Viscosity
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pT, // Temperature
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pF, // density, Rho
                                             std::shared_ptr<QuICC::DenseSM::IGenericProfile> pDF, // derivative of log(Rho)
                                             const QuICC::Equations::EquationParameters &eqParams, // physical nondimensional model parameters
                                             const MHDFloat c)
   {
      int nR = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int iR_;

      std::cerr << "nR = "<<nR<<" \n";

      auto nu        = pV->evaluateLP(r, 0, 0); 
      auto T         = pT->evaluateLP(r, 0, 0); 
      auto Rho       = pF->evaluateLP(r, 0, 0); 
      auto dLogRho   = pDF->evaluateLP(r, 0, 0);

      for(int iR = 0; iR < nR; ++iR)
      {
         iR_ = res.cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR);         
         int nTh = res.cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR); 

         std::cerr << "nTh = "<<nTh<<" \n";
         
         auto slice = computeViscousSlice(iR, iR_, c, v, Dv, nu(iR_), T(iR_), Rho(iR_), dLogRho(iR_));

         rS.addSlice(slice, iR);

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

         std::cerr <<" \n";
         std::cerr << "nPh = "<<nPh<<" \n";

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
         std::cerr << "v_r = "<<v.comp(FieldComponents::Physical::R).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_theta = "<<v.comp(FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_phi = "<<v.comp(FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
         std::cerr <<" \n";

         std::cerr << "v_rr = "<<Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::R).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_rtheta = "<<Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_rphi = "<<Dv.comp(FieldComponents::Physical::R, FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_thetar = "<<Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::R).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_thetatheta = "<<Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_thetaphi = "<<Dv.comp(FieldComponents::Physical::THETA, FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_phir = "<<Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::R).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_phitheta = "<<Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::THETA).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "v_phiphi = "<<Dv.comp(FieldComponents::Physical::PHI, FieldComponents::Physical::PHI).slice(iR).array()<<" \n";
         std::cerr <<" \n";
         std::cerr << "Q_nu = "<< (slice/c)*T(iR_)<<" \n";
         std::cerr <<" \n";
         std::cerr << "Di*Q_nu/T = "<< slice<<" \n";
         std::cerr <<" \n";
         std::cerr << "c = "<< c<<" \n";

         
      }

      
   // *** to print the result ** //
   // std::cerr << "NL(R) = "<<rS.data()<<" \n";

   
   }

}
}
