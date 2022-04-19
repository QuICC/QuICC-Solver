/**
 * @file Passthrough.cpp
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
#include "QuICC/PhysicalKernels/Passthrough.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

   Passthrough::Passthrough()
      : IPhysicalKernel()
   {
   }

   Passthrough::~Passthrough()
   {
   }

   std::size_t Passthrough::name() const
   {
      return this->mName;
   }

   void Passthrough::setField(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      IPhysicalKernel::setField(name, spField);
   }

   void Passthrough::setField(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      IPhysicalKernel::setField(name, spField);
   }

   void Passthrough::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      if(id == FieldComponents::Physical::SCALAR)
      {
         std::visit([&](auto&& p){rNLComp.setData(p->dom(0).phys().data());}, this->scalar(this->name()));
      } else
      {
         std::visit([&](auto&& p){rNLComp.setData(p->dom(0).phys().comp(id).data());}, this->vector(this->name()));
      }
   }

}
}
}
