/**
 * @file BenchmarkMagSteady.cpp
 * @brief Source of benchmark state for steady dynamo solution for magnetic field generator kernel
 */

// System includes
//

// Project includes
//
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkMagSteady.hpp"
#include "Types/Math.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   void BenchmarkMagSteady::init()
   {
   }

   void BenchmarkMagSteady::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
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
               amplitude = 4.0*std::sqrt(3.0)/Math::PI;
               funcR = (1.0 - (3.0/5.0)*std::pow(r,2));
               funcTh = std::cos(theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
            } else if(id == FieldComponents::Physical::THETA)
            {
               amplitude = -2.0*std::sqrt(3.0)/Math::PI;
               funcR = (2.0 - (12.0/5.0)*std::pow(r,2));
               funcTh = std::cos(theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPhCS,iTh,iR);
            } else if(id == FieldComponents::Physical::PHI)
            {
               amplitude = 20.0*std::sqrt(3.0)/Math::PI;
               funcR = std::pow(r,2)*(-1.0 + std::pow(r,2));
               funcTh = std::cos(theta)*std::sin(theta);
               rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
            }

         }
      }
   }

}
}
}
}
