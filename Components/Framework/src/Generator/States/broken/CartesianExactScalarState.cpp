/**
 * @file CartesianExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in cartesian geometries
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/CartesianExactScalarState.hpp"

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Math.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/NonDimensional/Lower1d.hpp"
#include "QuICC/NonDimensional/Upper1d.hpp"
#include "QuICC/NonDimensional/Lower2d.hpp"
#include "QuICC/NonDimensional/Upper2d.hpp"
#include "QuICC/NonDimensional/Lower3d.hpp"
#include "QuICC/NonDimensional/Upper3d.hpp"

namespace QuICC {

namespace Equations {

   CartesianExactScalarState::CartesianExactScalarState(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme)
      : IScalarEquation(spEqParams,spScheme), mTypeId(CartesianExactStateIds::CONSTANT), mModeA(3), mModeK(3)
   {
   }

   CartesianExactScalarState::~CartesianExactScalarState()
   {
   }

   void CartesianExactScalarState::setIdentity(const std::size_t name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void CartesianExactScalarState::setStateType(const CartesianExactStateIds::Id id)
   {
      this->mTypeId = id;
   }

   void CartesianExactScalarState::setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2)
   {
      this->mModeA.resize(2);
      this->mModeA(0) = a1;
      this->mModeA(1) = a2;

      this->mModeK.resize(2);
      this->mModeK(0) = k1;
      this->mModeK(1) = k2;
   }

   void CartesianExactScalarState::setModeOptions(const MHDFloat a1, const MHDFloat k1, const MHDFloat a2, const MHDFloat k2, const MHDFloat a3, const MHDFloat k3)
   {
      this->mModeA(0) = a1;
      this->mModeA(1) = a2;
      this->mModeA(2) = a3;

      this->mModeK(0) = k1;
      this->mModeK(1) = k2;
      this->mModeK(2) = k3;
   }

   void CartesianExactScalarState::setCoupling()
   {
      auto features = defaultCouplingFeature();
      features.at(CouplingFeature::Nonlinear) = true;
      features.at(CouplingFeature::AllowExplicit) = false;

      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, features);
   }

   void CartesianExactScalarState::initNLKernel(const bool force)
   {
      // Initialize if empty or forced
      if(force || !this->mspNLKernel)
      {
         // Assert on scalar component is used
         assert(compId == FieldComponents::Physical::SCALAR);

         if(this->mTypeId == CartesianExactStateIds::CONSTANT)
         {

         } else if(this->mTypeId == CartesianExactStateIds::PEYRET1DA)
         {

         } else if(this->mTypeId == CartesianExactStateIds::BXHELICOIDAL)
         {

         } else if(this->mTypeId == CartesianExactStateIds::BYHELICOIDAL)
         {

         } else if(this->mTypeId == CartesianExactStateIds::CONSTANTFIELD)
         {

         } else if(this->mTypeId == CartesianExactStateIds::NULLFIELD)
         {

         } else if(this->mTypeId == CartesianExactStateIds::PLANFORMSQUARES)
         {

         } else if(this->mTypeId == CartesianExactStateIds::TEST_DIFFUSION)
         {

         } else if(this->mTypeId == CartesianExactStateIds::TEST_BIDIFFUSION)
         {

         } else if(this->mTypeId == CartesianExactStateIds::TEST_BIDIFFUSION_SPLIT)
         {

         } else
         {
         }
      }
   }

//   void CartesianExactScalarState::computeNonlinear(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id compId) const
//   {
//      // Assert on scalar component is used
//      assert(compId == FieldComponents::Physical::SCALAR);
//
//      if(this->mTypeId == CartesianExactStateIds::CONSTANT)
//      {
//         rNLComp.rData().setConstant(this->mModeA.prod());
//
//      } else if(this->mTypeId == CartesianExactStateIds::PEYRET1DA)
//      {
//         #ifdef QUICC_SPATIALDIMENSION_3D
//            Array gI, gJ, gK;
//            this->buildGrid(gI, gJ, gK);
//         #else
//            Array gI, g;
//            this->buildGrid(gI, gJ);
//         #endif //QUICC_SPATIALDIMENSION_3D
//
//         MHDFloat k_;
//         //MHDFloat j_;
//         //MHDFloat i_;
//         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//         for(int iK = 0; iK < nK; ++iK)
//         {
//            k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//            for(int iJ = 0; iJ < nJ; ++iJ)
//            {
//               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//               //j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
//               for(int iI = 0; iI < nI; ++iI)
//               {
//                  //i_ = gI(iI);
//
//                  MHDFloat val = (-1.0 + 2.0*Math::PI)*std::cos(Math::PI*k_/2.0)/(2.0*Math::PI);
//
//                  rNLComp.setPoint(val, iI, iJ, iK);
//               }
//            }
//         }
//
//      } else if(this->mTypeId == CartesianExactStateIds::BXHELICOIDAL)
//      {
////Helicoidal Bx (Stellmach & Hansen, 2004)
//// written based on the PEYRET1DA case above
//// @ author: Stefano Maffei </maffei.ste@gmail.com/>
//// personal notes:
//// hope the presence of "this->" is not a problem (different from ShellExactVectorState.cpp)
//
//         Array gI, gJ, gK;
//         this->buildGrid(gI, gJ, gK);
//         MHDFloat k_;
//         //MHDFloat j_;
//         //MHDFloat i_;
//         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//         for(int iK = 0; iK < nK; ++iK)
//         {
//            k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//            for(int iJ = 0; iJ < nJ; ++iJ)
//            {
//               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//               //j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
//               for(int iI = 0; iI < nI; ++iI)
//               {
//                  //i_ = gI(iI);
//// If everything is what I think it is, I need only to modify this:
//                  MHDFloat val = -(std::sqrt(2)/2)*std::sin(Math::PI*k_/2.0) - (std::sqrt(2)/2)*std::sin(Math::PI*3.0*k_/2.0);
//
//                  rNLComp.setPoint(val, iI, iJ, iK);
//               }
//            }
//         }
//
//// ****** End of definition of BXHELICOIDAL ****** //
//      } else if(this->mTypeId == CartesianExactStateIds::BYHELICOIDAL)
//      {
////Helicoidal By (Stellmach & Hansen, 2004)
//
//         Array gI, gJ, gK;
//         this->buildGrid(gI, gJ, gK);
//         MHDFloat k_;
//         //MHDFloat j_;
//         //MHDFloat i_;
//         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//         for(int iK = 0; iK < nK; ++iK)
//         {
//            k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//            for(int iJ = 0; iJ < nJ; ++iJ)
//            {
//               //j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
//               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//               for(int iI = 0; iI < nI; ++iI)
//               {
//                  //i_ = gI(iI);
//                  MHDFloat val = std::sin(Math::PI*k_/2.0);
//                                    rNLComp.setPoint(val, iI, iJ, iK);
//               }
//            }
//         }
//
//
//// ****** End of definition of BYHELICOIDAL ****** //
//
//      } else if(this->mTypeId == CartesianExactStateIds::CONSTANTFIELD)
//      {
//// Simple constant (unity) field, for testing purposes
//         Array gI, gJ, gK;
//         this->buildGrid(gI, gJ, gK);
//         MHDFloat k_;
//
//
//         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//         for(int iK = 0; iK < nK; ++iK)
//         {
//            k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//            for(int iJ = 0; iJ < nJ; ++iJ)
//            {
//               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//               for(int iI = 0; iI < nI; ++iI)
//               {
//
//                  MHDFloat val = 1;
//
//                  rNLComp.setPoint(val, iI, iJ, iK);
//               }
//            }
//         }
//// ****** End of definition of a constant field Bx  ****** //
//
//      } else if(this->mTypeId == CartesianExactStateIds::NULLFIELD)
//      {
//// Simple constant null field, for testing purposes
//         Array gI, gJ, gK;
//         this->buildGrid(gI, gJ, gK);
//         MHDFloat k_;
//
//
//         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//         for(int iK = 0; iK < nK; ++iK)
//         {
//            k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//            for(int iJ = 0; iJ < nJ; ++iJ)
//            {
//               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//               for(int iI = 0; iI < nI; ++iI)
//               {
//
//                  MHDFloat val = 0;
//
//                  rNLComp.setPoint(val, iI, iJ, iK);
//               }
//            }
//         }
//
//// ****** End of definition of a null field  ****** //
//
//      } else if(this->mTypeId == CartesianExactStateIds::PLANFORMSQUARES)
//      {
//         #ifdef QUICC_SPATIALDIMENSION_3D
//            Array gI, gJ, gK;
//            this->buildGrid(gI, gJ, gK);
//         #else
//            int nK = 1;
//            int nJ = this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
//            int nI = this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
//
//            Array gK = -42*Array::Ones(1);
//            Array& gJ = this->mspMesh->at(0);
//            Array& gI = this->mspMesh->at(1);
//         #endif //QUICC_SPATIALDIMENSION_3D
//
//         int k_;
//         int j_;
//         int i_;
//         MHDFloat gk_;
//         MHDFloat gj_;
//         MHDFloat gi_;
//         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//         for(int iK = 0; iK < nK; ++iK)
//         {
//            k_ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK);
//            gk_ = gK(k_);
//            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//            for(int iJ = 0; iJ < nJ; ++iJ)
//            {
//               j_ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK);
//               gj_ = gJ(j_);
//               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//               for(int iI = 0; iI < nI; ++iI)
//               {
//                  i_ = iI;
//                  gi_ = gI(i_);
//
//                  MHDFloat val = (1+gk_)*(1.0 + std::cos(this->mModeK(0)*gj_) + std::cos(this->mModeK(1)*gi_));
//
//                  rNLComp.setPoint(val, iI, iJ, iK);
//               }
//            }
//         }
//
//      } else if(this->mTypeId == CartesianExactStateIds::TEST_DIFFUSION)
//      {
//         #if defined QUICC_SPATIALDIMENSION_3D
//            Array gI, gK, gJ;
//            this->buildGrid(gI, gJ, gK);
//            #if defined QUICC_SPATIALSCHEME_FFF
//            #elif defined QUICC_SPATIALSCHEME_TFF
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper1d::id());
//            #elif defined QUICC_SPATIALSCHEME_TFT
//               MHDFloat xi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat xo = this->eqParams().nd(NonDimensional::Upper1d::id());
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower3d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper3d::id());
//            #elif defined QUICC_SPATIALSCHEME_TTT
//               MHDFloat xi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat xo = this->eqParams().nd(NonDimensional::Upper1d::id());
//               MHDFloat yi = this->eqParams().nd(NonDimensional::Lower2d::id());
//               MHDFloat yo = this->eqParams().nd(NonDimensional::Upper2d::id());
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower3d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper3d::id());
//            #endif //defined QUICC_SPATIALSCHEME_FFF
//
//            MHDFloat k_;
//            MHDFloat j_;
//            MHDFloat i_;
//            int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//            for(int iK = 0; iK < nK; ++iK)
//            {
//               k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//               int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//               for(int iJ = 0; iJ < nJ; ++iJ)
//               {
//                  j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
//                  int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//                  for(int iI = 0; iI < nI; ++iI)
//                  {
//                     i_ = gI(iI);
//
//                     MHDFloat PI = Math::PI;
//                     MHDFloat funcT = (-1.0 + 2.0*Math::PI)/(2.0*Math::PI);
//
//                     #if defined QUICC_SPATIALSCHEME_FFF
//                        MHDFloat val = 0.0;
//                        for(int kx = 0; kx < 3; ++kx)
//                        {
//                           for(int ky = 0; ky < 3; ++ky)
//                           {
//                              for(int kz = 0; kz < 3; ++kz)
//                              {
//                                 MHDFloat funcX = (std::cos(kx*k_) + std::sin(kx*k_));
//                                 MHDFloat funcY = (std::cos(ky*j_) + std::sin(ky*j_));
//                                 MHDFloat funcZ = (std::cos(kz*i_) + std::sin(kz*i_));
//                                 val += funcT*funcX*funcY*funcZ;
//                              }
//                           }
//                        }
//                     #elif defined QUICC_SPATIALSCHEME_TFF
//                        MHDFloat val = 0.0;
//                        for(int kx = 0; kx < 3; ++kx)
//                        {
//                           for(int ky = 0; ky < 3; ++ky)
//                           {
//                              MHDFloat funcX = (std::cos(kx*j_) + std::sin(kx*j_));
//                              MHDFloat funcY = (std::cos(ky*i_) + std::sin(ky*i_));
//                              MHDFloat z = (2.0*k_ - (zo + zi))/(zo - zi);
//                              MHDFloat funcZ = (1.0 - std::pow(z,2))*std::sin(PI*z);
//                              val += funcT*funcX*funcY*funcZ;
//                           }
//                        }
//                     #elif defined QUICC_SPATIALSCHEME_TFT
//                        MHDFloat x = (2.0*k_ - (xo + xi))/(xo - xi);
//                        MHDFloat z = (2.0*i_ - (zo + zi))/(zo - zi);
//                        MHDFloat val = 0.0;
//                        for(int ky = 0; ky < 3; ++ky)
//                        {
//                           MHDFloat funcX = (1.0 - std::pow(x,2))*std::sin(PI*x);
//                           MHDFloat funcY = (std::cos(ky*j_) + std::sin(ky*j_));
//                           MHDFloat funcZ = (1.0 - std::pow(z,2))*std::sin(PI*z);
//                           val += funcT*funcX*funcY*funcZ;
//                        }
//                     #elif defined QUICC_SPATIALSCHEME_TTT
//                        MHDFloat x = (2.0*k_ - (xo + xi))/(xo - xi);
//                        MHDFloat funcX = (1.0 - std::pow(x,2))*std::sin(PI*x);
//                        MHDFloat y = (2.0*j_ - (yo + yi))/(yo - yi);
//                        MHDFloat funcY = (1.0 - std::pow(y,2))*std::sin(PI*y);
//                        MHDFloat z = (2.0*i_ - (zo + zi))/(zo - zi);
//                        MHDFloat funcZ = (1.0 - std::pow(z,2))*std::sin(PI*z);
//                        MHDFloat val = funcT*funcX*funcY*funcZ;
//                     #endif //defined QUICC_SPATIALSCHEME_FFF
//
//                     rNLComp.setPoint(val, iI, iJ, iK);
//                  }
//               }
//            }
//
//         #elif defined QUICC_SPATIALDIMENSION_2D
//            // NOT IMPLEMENTED YET
//         #elif defined QUICC_SPATIALDIMENSION_1D
//            // NOT IMPLEMENTED YET
//         #endif //defined QUICC_SPATIALDIMENSION_3D
//
//      } else if(this->mTypeId == CartesianExactStateIds::TEST_BIDIFFUSION)
//      {
//         #if defined QUICC_SPATIALDIMENSION_3D
//            Array gI, gK, gJ;
//            this->buildGrid(gI, gJ, gK);
//            #if defined QUICC_SPATIALSCHEME_FFF
//            #elif defined QUICC_SPATIALSCHEME_TFF
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper1d::id());
//            #elif defined QUICC_SPATIALSCHEME_TFT
//               MHDFloat xi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat xo = this->eqParams().nd(NonDimensional::Upper1d::id());
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower3d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper3d::id());
//            #elif defined QUICC_SPATIALSCHEME_TTT
//               MHDFloat xi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat xo = this->eqParams().nd(NonDimensional::Upper1d::id());
//               MHDFloat yi = this->eqParams().nd(NonDimensional::Lower2d::id());
//               MHDFloat yo = this->eqParams().nd(NonDimensional::Upper2d::id());
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower3d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper3d::id());
//            #endif //defined QUICC_SPATIALSCHEME_FFF
//
//            MHDFloat k_;
//            MHDFloat j_;
//            MHDFloat i_;
//            int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//            for(int iK = 0; iK < nK; ++iK)
//            {
//               k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//               int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//               for(int iJ = 0; iJ < nJ; ++iJ)
//               {
//                  j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
//                  int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//                  for(int iI = 0; iI < nI; ++iI)
//                  {
//                     i_ = gI(iI);
//
//                     MHDFloat PI = Math::PI;
//                     MHDFloat funcT = (-1.0 + 2.0*Math::PI)/(2.0*Math::PI);
//
//                     #if defined QUICC_SPATIALSCHEME_FFF
//                        MHDFloat val = 0.0;
//                        for(int kx = 0; kx < 3; ++kx)
//                        {
//                           for(int ky = 0; ky < 3; ++ky)
//                           {
//                              for(int kz = 0; kz < 3; ++kz)
//                              {
//                                 MHDFloat funcX = (std::cos(kx*k_) + std::sin(kx*k_));
//                                 MHDFloat funcY = (std::cos(ky*j_) + std::sin(ky*j_));
//                                 MHDFloat funcZ = (std::cos(kz*i_) + std::sin(kz*i_));
//                                 val += funcT*funcX*funcY*funcZ;
//                              }
//                           }
//                        }
//                     #elif defined QUICC_SPATIALSCHEME_TFF
//                        MHDFloat val = 0.0;
//                        for(int kx = 0; kx < 3; ++kx)
//                        {
//                           for(int ky = 0; ky < 3; ++ky)
//                           {
//                              MHDFloat funcX = (std::cos(kx*j_) + std::sin(kx*j_));
//                              MHDFloat funcY = (std::cos(ky*i_) + std::sin(ky*i_));
//                              MHDFloat z = (2.0*k_ - (zo + zi))/(zo - zi);
//                              MHDFloat funcZ = std::pow(1.0 - std::pow(z,2),2)*std::sin(PI*z);
//                              val += funcT*funcX*funcY*funcZ;
//                           }
//                        }
//                     #elif defined QUICC_SPATIALSCHEME_TFT
//                        MHDFloat x = (2.0*k_ - (xo + xi))/(xo - xi);
//                        MHDFloat z = (2.0*i_ - (zo + zi))/(zo - zi);
//                        MHDFloat val = 0.0;
//                        for(int ky = 0; ky < 3; ++ky)
//                        {
//                           MHDFloat funcX = (1.0 - std::pow(x,2))*std::sin(PI*x);
//                           MHDFloat funcY = (std::cos(ky*j_) + std::sin(ky*j_));
//                           MHDFloat funcZ = (1.0 - std::pow(z,2))*std::sin(PI*z);
//                           val += funcT*funcX*funcY*funcZ;
//                        }
//                     #elif defined QUICC_SPATIALSCHEME_TTT
//                        MHDFloat x = (2.0*k_ - (xo + xi))/(xo - xi);
//                        MHDFloat funcX = (1.0 - std::pow(x,2))*std::sin(PI*x);
//                        MHDFloat y = (2.0*j_ - (yo + yi))/(yo - yi);
//                        MHDFloat funcY = (1.0 - std::pow(y,2))*std::sin(PI*y);
//                        MHDFloat z = (2.0*i_ - (zo + zi))/(zo - zi);
//                        MHDFloat funcZ = (1.0 - std::pow(z,2))*std::sin(PI*z);
//                        MHDFloat val = funcT*funcX*funcY*funcZ;
//                     #endif //defined QUICC_SPATIALSCHEME_FFF
//
//                     rNLComp.setPoint(val, iI, iJ, iK);
//                  }
//               }
//            }
//
//         #elif defined QUICC_SPATIALDIMENSION_2D
//            // NOT IMPLEMENTED YET
//         #elif defined QUICC_SPATIALDIMENSION_1D
//            // NOT IMPLEMENTED YET
//         #endif //defined QUICC_SPATIALDIMENSION_3D
//
//      } else if(this->mTypeId == CartesianExactStateIds::TEST_BIDIFFUSION_SPLIT)
//      {
//         MHDFloat PI = Math::PI;
//         #if defined QUICC_SPATIALDIMENSION_3D
//            Array gI, gK, gJ;
//            this->buildGrid(gI, gJ, gK);
//            #if defined QUICC_SPATIALSCHEME_FFF
//            #elif defined QUICC_SPATIALSCHEME_TFF
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper1d::id());
//               MHDFloat dz = zo - zi;
//            #elif defined QUICC_SPATIALSCHEME_TFT
//               MHDFloat xi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat xo = this->eqParams().nd(NonDimensional::Upper1d::id());
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower3d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper3d::id());
//            #elif defined QUICC_SPATIALSCHEME_TTT
//               MHDFloat xi = this->eqParams().nd(NonDimensional::Lower1d::id());
//               MHDFloat xo = this->eqParams().nd(NonDimensional::Upper1d::id());
//               MHDFloat yi = this->eqParams().nd(NonDimensional::Lower2d::id());
//               MHDFloat yo = this->eqParams().nd(NonDimensional::Upper2d::id());
//               MHDFloat zi = this->eqParams().nd(NonDimensional::Lower3d::id());
//               MHDFloat zo = this->eqParams().nd(NonDimensional::Upper3d::id());
//            #endif //defined QUICC_SPATIALSCHEME_FFF
//
//            MHDFloat k_;
//            MHDFloat j_;
//            MHDFloat i_;
//            int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//            for(int iK = 0; iK < nK; ++iK)
//            {
//               k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//               int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//               for(int iJ = 0; iJ < nJ; ++iJ)
//               {
//                  j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
//                  int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//                  for(int iI = 0; iI < nI; ++iI)
//                  {
//                     i_ = gI(iI);
//
//
//                     #if defined QUICC_SPATIALSCHEME_FFF
//                     #elif defined QUICC_SPATIALSCHEME_TFF
//                        MHDFloat val = 0.0;
//                        MHDFloat funcT = (-1.0 + 2.0*PI)/(2.0*Math::PI*std::pow(dz,2));
//                        MHDFloat z = (2.0*k_ - (zo + zi))/(zo - zi);
//                        MHDFloat cz = std::cos(PI*z);
//                        MHDFloat sz = std::sin(PI*z);
//                        MHDFloat zz = (-1.0+std::pow(z,2));
//                        for(int kx = 0; kx < 3; ++kx)
//                        {
//                           MHDFloat funcX = (std::cos(kx*j_) + std::sin(kx*j_));
//                           for(int ky = 0; ky < 3; ++ky)
//                           {
//                              MHDFloat kk = kx*kx + ky*ky;
//                              MHDFloat funcY = (std::cos(ky*i_) + std::sin(ky*i_));
//                              MHDFloat funcZ = (32.0*PI*z*zz*cz-(16.0-48.0*std::pow(z,2)+std::pow(dz,2)*kk*std::pow(zz,2)+4.0*std::pow(PI,2)*std::pow(zz,2))*sz);
//                              val += funcT*funcX*funcY*funcZ;
//                           }
//                        }
//                     #elif defined QUICC_SPATIALSCHEME_TFT
//                     #elif defined QUICC_SPATIALSCHEME_TTT
//                     #endif //defined QUICC_SPATIALSCHEME_FFF
//
//                     rNLComp.setPoint(val, iI, iJ, iK);
//                  }
//               }
//            }
//
//         #elif defined QUICC_SPATIALDIMENSION_2D
//            // NOT IMPLEMENTED YET
//         #elif defined QUICC_SPATIALDIMENSION_1D
//            // NOT IMPLEMENTED YET
//         #endif //defined QUICC_SPATIALDIMENSION_3D
//
//      } else
//      {
//         #ifdef QUICC_SPATIALDIMENSION_3D
//            Array gI, gJ, gK;
//            this->buildGrid(gI, gJ, gK);
//         #else
//            int nK = 1;
//            int nJ = this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
//            int nI = this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
//
//            Array gK = -42*Array::Ones(1);
//            Array gJ = this->mspMesh->at(0);
//            Array gI = this->mspMesh->at(1);
//         #endif //QUICC_SPATIALDIMENSION_3D
//
//         MHDFloat k_;
//         MHDFloat j_;
//         MHDFloat i_;
//         int nK = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT3D>();
//         for(int iK = 0; iK < nK; ++iK)
//         {
//            k_ = gK(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT3D>(iK));
//            int nJ = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DAT2D>(iK);
//            for(int iJ = 0; iJ < nJ; ++iJ)
//            {
//               int nI = this->res().cpu()->dim(Dimensions::Transform::TRAND)->dim<Dimensions::Data::DATF1D>(iK);
//               j_ = gJ(this->res().cpu()->dim(Dimensions::Transform::TRAND)->idx<Dimensions::Data::DAT2D>(iJ, iK));
//               for(int iI = 0; iI < nI; ++iI)
//               {
//                  i_ = gI(iI);
//
//                  MHDFloat val = 0.0;
//                  if(static_cast<int>(this->mTypeId) < 100)
//                  {
//                     Array grid(3);
//                     grid(0) = k_;
//                     grid(1) = j_;
//                     grid(2) = i_;
//                     val = CartesianExactStateIds::exact3D(this->mTypeId, this->mModeA, this->mModeK, grid);
//                  } else if(static_cast<int>(this->mTypeId) >= 100)
//                  {
//                     Array grid(2);
//                     grid(0) = j_;
//                     grid(1) = i_;
//                     val = CartesianExactStateIds::exact2D(this->mTypeId, this->mModeA, this->mModeK, grid);
//                  }
//
//                  rNLComp.setPoint(val, iI, iJ, iK);
//               }
//            }
//         }
//      }
//   }

   MHDVariant CartesianExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
   {
     // Assert on scalar component is used
     assert(compId == FieldComponents::Spectral::SCALAR);

     return MHDComplex(0);
   }

   void CartesianExactScalarState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::After::id());

      // Get reference to spatial scheme
      const auto& ss = this->ss();

      // Add unknown to requirements: is scalar?
      auto& req = this->mRequirements.addField(this->name(), FieldRequirement(true, ss.spectral(), ss.physical()));
      req.enableSpectral();
      req.enablePhysical();
   }

   void CartesianExactScalarState::buildGrid(Array& g1D, Array& g2D, Array& g3D) const
   {
      int nK = this->res().sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ = this->res().sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nI = this->res().sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      g3D = this->mspMesh->at(2);
      g2D = this->mspMesh->at(1);
      g1D = this->mspMesh->at(0);
   }

}
}
