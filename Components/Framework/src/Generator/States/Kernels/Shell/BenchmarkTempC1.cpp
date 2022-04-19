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
#include "QuICC/Generator/States/Kernels/Shell/BenchmarkTempC1.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Shell {

   BenchmarkTempC1::BenchmarkTempC1()
      : IPhysicalKernel()
   {
   }

   BenchmarkTempC1::~BenchmarkTempC1()
   {
   }

   void BenchmarkTempC1::init(const MHDFloat ri, const MHDFloat ro)
   {
      this->mRi = ri;
      this->mRo = ro;
   }

   void BenchmarkTempC1::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);

      MHDFloat ri = this->mRi;
      MHDFloat ro = this->mRo;
      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);

      MHDFloat funcR;
      MHDFloat funcTh;
      Array funcPh = (4.0*phGrid).array().cos();

      MHDFloat amplitude = 21.0/std::sqrt(17920*Math::PI);

      MHDFloat r;
      MHDFloat x;
      MHDFloat theta;

      rNLComp.rData().setConstant(0);
      nR = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
         x = 2.0*r - ri - ro;
         funcR = 1.0 - 3.0*std::pow(x,2) + 3.0*std::pow(x,4) - std::pow(x,6);

         nTh = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
            funcTh = std::pow(std::sin(theta),4);

            rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
         }
      }
   }

}
}
}
}
