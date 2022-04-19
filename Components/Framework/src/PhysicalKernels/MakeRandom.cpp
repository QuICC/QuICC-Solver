/**
 * @file MakeRandom.cpp
 * @brief Source of trivial kernel to make field random with given scale
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalKernels/MakeRandom.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

   MakeRandom::MakeRandom()
      : IPhysicalKernel()
   {
   }

   MakeRandom::~MakeRandom()
   {
   }

   void MakeRandom::init(const MHDFloat scale)
   {
      this->mScale = scale;
   }

   void MakeRandom::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      rNLComp.rData().setRandom();
      rNLComp.rData() *= this->mScale;
   }

}
}
}
