/**
 * @file BenchmarkMagC1.cpp
 * @brief Source of benchmark state C1 for magnetic field generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Shell/BenchmarkMagC1.hpp"

// Project includes
//
#include "Types/Constants.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Shell {

   BenchmarkMagC1::BenchmarkMagC1()
      : IPhysicalKernel()
   {
   }

   BenchmarkMagC1::~BenchmarkMagC1()
   {
   }

   void BenchmarkMagC1::init(const MHDFloat ri, const MHDFloat ro)
   {
      this->mRi = ri;
      this->mRo = ro;
   }

   void BenchmarkMagC1::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      MHDFloat ri = this->mRi;
      MHDFloat ro = this->mRo;
      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);

      Array funcPh = Array::Ones(nPh);
      MHDFloat funcR = 1.0;
      MHDFloat funcTh = 1.0;
      MHDFloat amplitude = 1.0;

      MHDFloat r;
      MHDFloat theta;
      MHDFloat scale = 1.0/std::sqrt(2.0);
      nR = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
         nTh = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            if(id == FieldComponents::Physical::R)
            {
               amplitude = scale*5.0/8.0;
               funcR = 8.0*ro - 6.0*r - 2.0*std::pow(ri,4)/std::pow(r,3);
               funcTh = std::cos(theta);
            } else if(id == FieldComponents::Physical::THETA)
            {
               amplitude = -scale*5.0/8.0;
               funcR = 8.0*ro - 9.0*r + std::pow(ri,4)/std::pow(r,3);
               funcTh = std::sin(theta);
            } else if(id == FieldComponents::Physical::PHI)
            {
               amplitude = scale*5.0;
               funcR = std::sin(Math::PI*(r - ri));
               funcTh = std::sin(2*theta);
            }

            rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
         }
      }
   }

}
}
}
}
