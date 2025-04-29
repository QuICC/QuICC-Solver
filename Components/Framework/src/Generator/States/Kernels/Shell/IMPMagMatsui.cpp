/**
 * @file IMPMagMatsui.cpp
 * @brief Source of IMP Magnetic state from Matsui generator kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/States/Kernels/Shell/IMPMagMatsui.hpp"

// Project includes
//
#include "Types/Math.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

namespace Shell {

   IMPMagMatsui::IMPMagMatsui()
      : IPhysicalKernel()
   {
   }

   IMPMagMatsui::~IMPMagMatsui()
   {
   }

   void IMPMagMatsui::init(const MHDFloat ri, const MHDFloat ro)
   {
      this->mRi = ri;
      this->mRo = ro;
   }

   void IMPMagMatsui::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Initialize to zero
      rNLComp.rData().setZero();

      int nR = this->spRes()->sim().dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nTh = this->spRes()->sim().dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nPh = this->spRes()->sim().dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      //MHDFloat ri = this->mRi;
      //MHDFloat ro = this->mRo;
      //Array& rGrid = this->mspMesh->at(0);
      Array& thGrid = this->mspMesh->at(1);
      //Array& phGrid = this->mspMesh->at(2);

      Array funcPh = Array::Ones(nPh);
      MHDFloat funcR = 1.0;
      MHDFloat funcTh = 1.0;
      MHDFloat amplitude = 1.0;

      MHDFloat theta;
      nR = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iR = 0; iR < nR; ++iR)
      {
         //MHDFloat r = rGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
         nTh = this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            theta = thGrid(this->spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

            if(id == FieldComponents::Physical::R)
            {
               amplitude = 1.0;
               funcR = 1.0;
               funcTh = std::cos(theta);
            } else if(id == FieldComponents::Physical::THETA)
            {
               amplitude = -1.0;
               funcR = 1.0;
               funcTh = std::sin(theta);
            } else if(id == FieldComponents::Physical::PHI)
            {
               amplitude = 0.0;
               funcR = 0.0;
               funcTh = 0.0;
            }

            rNLComp.addProfile(amplitude*funcR*funcTh*funcPh,iTh,iR);
         }
      }
   }

}
}
}
}
