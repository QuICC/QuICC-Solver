/**
 * @file MakeConstant.cpp
 * @brief Source of trivial kernel to make field constant
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalKernels/MakeConstant.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

   MakeConstant::MakeConstant()
      : IPhysicalKernel()
   {
   }

   MakeConstant::~MakeConstant()
   {
   }

   void MakeConstant::init(const MHDFloat value)
   {
      this->mValue = value;
   }

   void MakeConstant::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      rNLComp.rData().setConstant(this->mValue);
   }

}
}
}
