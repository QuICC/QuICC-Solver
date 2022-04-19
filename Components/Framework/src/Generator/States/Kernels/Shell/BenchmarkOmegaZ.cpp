/**
 * @file BenchmarkOmegaZ.cpp
 * @brief Source of benchmark Omega Z velocity field generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Shell/BenchmarkOmegaZ.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Shell {

   BenchmarkOmegaZ::BenchmarkOmegaZ()
      : IPhysicalKernel()
   {
   }

   BenchmarkOmegaZ::~BenchmarkOmegaZ()
   {
   }

   void BenchmarkOmegaZ::init(const MHDFloat ri, const MHDFloat ro)
   {
      this->mRi = ri;
      this->mRo = ro;
   }

   void BenchmarkOmegaZ::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);

      Array funcPh = Array::Ones(nPh);
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
               funcR = 0.0;
               funcTh = 0.0;
            } else if(id == FieldComponents::Physical::THETA)
            {
               amplitude = 0.0;
               funcR = 0.0;
               funcTh = 0.0;
            } else if(id == FieldComponents::Physical::PHI)
            {
               amplitude = 1.0;
               funcR = r;
               funcTh = std::sin(theta);
            }

            rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
    	   }
    	}
   }

}
}
}
}
