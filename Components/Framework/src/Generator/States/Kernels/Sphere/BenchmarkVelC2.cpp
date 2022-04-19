/**
 * @file BenchmarkVelC2.cpp
 * @brief Source of benchmark state C2 for velocity field generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkVelC2.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   BenchmarkVelC2::BenchmarkVelC2()
      : IPhysicalKernel()
   {
   }

   BenchmarkVelC2::~BenchmarkVelC2()
   {
   }

   void BenchmarkVelC2::init()
   {
   }

   void BenchmarkVelC2::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      Array& phGrid = this->mspMesh->at(2);

      Array funcPh = Array::Ones(nPh);
      Array funcPhC1 = (phGrid).array().cos();
      Array funcPhS1 = (phGrid).array().sin();
      MHDFloat funcR = 1.0;
      MHDFloat amplitude = 1.0;
      MHDFloat funcTh = 1.0;

      MHDFloat r;
      MHDFloat theta;
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
               amplitude = 0.0;
               rNLComp.addProfile(amplitude*funcPh,iTh,iR);
            } else if(id == FieldComponents::Physical::THETA)
            {
               amplitude = -10.0/(7.0*std::sqrt(3.0));
               funcTh = std::cos(theta);

               funcR = 3.0*std::pow(r,2)*(-147.0 + 343.0*std::pow(r,2) - 217.0*std::pow(r,4) + 29.0*std::pow(r,6));
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPhC1,iTh,iR);

               funcR = 14.0*std::pow(r,2)*(-9.0 - 125.0*std::pow(r,2) + 39.0*std::pow(r,4) + 27.0*std::pow(r,6));
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPhS1,iTh,iR);
            } else if(id == FieldComponents::Physical::PHI)
            {
               amplitude = -5.0/5544.;

               funcR = 7.0*r*(43700.0 - 58113.0*std::pow(r,2) - 15345.0*std::pow(r,4) + 1881.0*std::pow(r,6) + 20790.0*std::pow(r,8));
               funcTh = std::sin(theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);

               funcR = 7.0*1485*std::pow(r,3)*(-9.0 + 115.0*std::pow(r,2) - 167.0*std::pow(r,4) + 70.0*std::pow(r,6));
               funcTh = std::sin(3.0*theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);

               funcR = 528*std::sqrt(3)*std::pow(r,2)*14*(-9.0 - 125.0*std::pow(r,2) + 39.0*std::pow(r,4) + 27.0*std::pow(r,6));
               funcTh = std::cos(2.0*theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPhC1,iTh,iR);

               funcR = 528*std::sqrt(3)*std::pow(r,2)*3*(147.0 - 343.0*std::pow(r,2) + 217.0*std::pow(r,4) - 29.0*std::pow(r,6));
               funcTh = std::cos(2.0*theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPhS1,iTh,iR);
            }
         }
      }

   }

}
}
}
}
