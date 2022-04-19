/**
 * @file DoNothing.cpp
 * @brief Source of trivial passthrough kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalKernels/DoNothing.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

   DoNothing::DoNothing()
      : IPhysicalKernel()
   {
   }

   DoNothing::~DoNothing()
   {
   }

   std::size_t DoNothing::name() const
   {
      return this->mName;
   }

   void DoNothing::setField(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      IPhysicalKernel::setField(name, spField);
   }

   void DoNothing::setField(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      IPhysicalKernel::setField(name, spField);
   }

   void DoNothing::compute(Framework::Selector::PhysicalScalarField&, FieldComponents::Physical::Id) const
   {
   }

}
}
}
