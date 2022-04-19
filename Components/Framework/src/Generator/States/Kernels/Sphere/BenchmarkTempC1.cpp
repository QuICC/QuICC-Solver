/**
 * @file BenchmarkTempC1.cpp
 * @brief Source of benchmark state C1 for temperature generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkTempC1.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   BenchmarkTempC1::BenchmarkTempC1()
      : IPhysicalKernel()
   {
   }

   BenchmarkTempC1::~BenchmarkTempC1()
   {
   }

   void BenchmarkTempC1::init(const MHDFloat amplitude_bg, const MHDFloat epsilon)
   {
      this->mBg = amplitude_bg;
      this->mEpsilon = epsilon;
   }

   void BenchmarkTempC1::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
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
      MHDFloat funcR3;
      MHDFloat funcTh0;
      MHDFloat funcTh3;
      Array funcPh0 = Array::Ones(nPh);
      Array funcPh3 = (3.0*phGrid).array().cos() + (3.0*phGrid).array().sin();

      // Background state is not solved for (no source term but background is imposed explicitly)
      MHDFloat amp0 = this->mBg;
      MHDFloat eps = this->mEpsilon;
      MHDFloat amp3 = (eps/8.0)*std::sqrt(35.0/Math::PI);

      MHDFloat r;
      MHDFloat theta;

      rNLComp.rData().setConstant(0);
      nR = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
         funcR0 = (1.0 - std::pow(r,2));
         funcR3 = std::pow(r,3)*(1.0 - std::pow(r,2));

         nTh = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
            funcTh0 = 1.0;
            funcTh3 = std::pow(std::sin(theta),3);

            rNLComp.addProfile(amp0*funcR0*funcTh0*funcPh0,iTh,iR);
            rNLComp.addProfile(amp3*funcR3*funcTh3*funcPh3,iTh,iR);
         }
      }
   }

}
}
}
}
