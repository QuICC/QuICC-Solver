/**
 * @file BenchmarkTempSteady.cpp
 * @brief Source of benchmark state steady dynamo solution for temperature generator kernel
 */

// System includes
//

// Project includes
//
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkTempSteady.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   void BenchmarkTempSteady::init(const MHDFloat amplitude_bg, const MHDFloat epsilon)
   {
      this->mBg = amplitude_bg;
      this->mEpsilon = epsilon;
   }

   void BenchmarkTempSteady::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);

      MHDFloat funcR0;
      MHDFloat funcR2;
      MHDFloat funcTh0;
      MHDFloat funcTh2;
      Array funcPh0 = Array::Ones(nPh);
      Array funcPh2 = (2.0*phGrid).array().cos();

      // Background state is not solved for (no source term but background is imposed explicitly)
      MHDFloat eps = (2*std::sqrt(10)/(5.0*Math::PI))*this->mEpsilon;

      MHDFloat r;
      MHDFloat theta;

      auto& tRes = *this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D);

      rNLComp.rData().setConstant(0);
      nR = tRes.dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         r = rGrid(tRes.idx<Dimensions::Data::DAT3D>(iR));
         funcR0 = (1.0 - std::pow(r,2));
         funcR2 = std::pow(r,2)*(1.0 - std::pow(r,2));

         nTh = tRes.dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(tRes.idx<Dimensions::Data::DAT2D>(iTh, iR));
            funcTh0 = 1.0;
            funcTh2 = std::pow(std::sin(theta),2);

            // background state
            rNLComp.addProfile(this->mBg*funcR0*funcTh0*funcPh0,iTh,iR);

            // perturbation
            rNLComp.addProfile(eps*funcR2*funcTh2*funcPh2,iTh,iR);
         }
      }
   }

} // Sphere
} // Kernel
} // Physical
} // QuICC
