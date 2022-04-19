/**
 * @file BenchmarkMagC2.cpp
 * @brief Source of benchmark state C2 for magnetic field generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkMagC2.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   BenchmarkMagC2::BenchmarkMagC2()
      : IPhysicalKernel()
   {
   }

   BenchmarkMagC2::~BenchmarkMagC2()
   {
   }

   void BenchmarkMagC2::init()
   {
   }

   void BenchmarkMagC2::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
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
      Array funcPhCS = (phGrid).array().cos() + (phGrid).array().sin();
      Array funcPhC_S = (phGrid).array().cos() - (phGrid).array().sin();
      MHDFloat funcR = 1.0;
      MHDFloat funcTh = 1.0;
      MHDFloat amplitude = 1.0;

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
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
            } else if(id == FieldComponents::Physical::THETA)
            {
               amplitude = -3.0/2.0;
               funcR = r*(-1.0 + 4.0*std::pow(r,2) - 6.0*std::pow(r,4) + 3.0*std::pow(r,6));
               funcTh = 1.0;
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPhCS,iTh,iR);
            } else if(id == FieldComponents::Physical::PHI)
            {
               amplitude = -3.0/4.0;

               funcR = 3.0*std::pow(r,2)*(-1.0 + std::pow(r,2))*(2.0 - 5.0*std::pow(r,2) + 4.0*std::pow(r,4));
               funcTh = std::cos(theta)*std::sin(theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);

               funcR = 2.0*r*(-1.0 + std::pow(r,2))*(1.0 - 3.0*std::pow(r,2) + 3.0*std::pow(r,4));
               funcTh = std::cos(theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPhC_S,iTh,iR);
            }

         }
      }
   }

}
}
}
}
