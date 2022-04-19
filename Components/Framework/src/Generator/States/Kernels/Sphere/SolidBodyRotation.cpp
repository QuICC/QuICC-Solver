/**
 * @file SolidBodyRotation.cpp
 * @brief Source of solid body rotation generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Sphere/SolidBodyRotation.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   SolidBodyRotation::SolidBodyRotation()
      : IPhysicalKernel(), mX(0.0), mY(0.0), mZ(0.0)
   {
   }

   SolidBodyRotation::~SolidBodyRotation()
   {
   }

   void SolidBodyRotation::init(const MHDFloat x, const MHDFloat y, const MHDFloat z)
   {
      this->mX = x;
      this->mY = y;
      this->mZ = z;
   }

   void SolidBodyRotation::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
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
               rNLComp.addProfile(0.0*funcPh,iTh,iR);
            } else if(id == FieldComponents::Physical::THETA)
            {
               funcTh = 1.0;
               funcR = r;

               rNLComp.addProfile(this->mY*funcR*funcTh*funcPhC1,iTh,iR);
               rNLComp.subProfile(this->mX*funcR*funcTh*funcPhS1,iTh,iR);
            } else if(id == FieldComponents::Physical::PHI)
            {
               funcR = r;

               funcTh = std::sin(theta);
               rNLComp.addProfile(this->mZ*funcR*funcTh*funcPh,iTh,iR);

               funcTh = std::cos(theta);
               rNLComp.subProfile(this->mX*funcR*funcTh*funcPhC1,iTh,iR);
               rNLComp.subProfile(this->mY*funcR*funcTh*funcPhS1,iTh,iR);
            }
         }
      }

   }

}
}
}
}
