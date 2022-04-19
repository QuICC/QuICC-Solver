/**
 * @file MalkusField.cpp
 * @brief Source of uniform axial field generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Sphere/MalkusField.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Sphere {

   MalkusField::MalkusField()
      : IPhysicalKernel(), mA(1.0)
   {
   }

   MalkusField::~MalkusField()
   {
   }

   void MalkusField::init(const MHDFloat amplitude)
   {
      this->mA = amplitude;
   }

   void MalkusField::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      //Array& phGrid = this->mspMesh->at(2);

      Array funcPh = Array::Ones(nPh);
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

            if(id == FieldComponents::Physical::PHI)
            {
               funcTh = std::sin(theta);
               rNLComp.addProfile(this->mA*r*funcTh*funcPh,iTh,iR);
            }

         }
      }
   }

}
}
}
}
