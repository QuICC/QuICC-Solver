/**
 * @file CurlKernel.cpp
 * @brief Source of physical space inertia kernel
 */

// System includes
//

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/CurlKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   std::size_t CurlKernel::name() const
   {
      return this->mName;
   }

   void CurlKernel::setVelocity(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      this->setField(name, spField);
   }

   void CurlKernel::init(const MHDFloat inertia)
   {
      // Set scaling constants
      this->mInertia = inertia;
   }

   void CurlKernel::setMesh(std::shared_ptr<std::vector<Array> > spMesh)
   {
      IPhysicalKernel::setMesh(spMesh);

      this->mRadius = this->mspMesh->at(0);
   }

   void CurlKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      std::visit([&](auto&& v)
      {
         rNLComp.setData(v->dom(0).curl().comp(id).data());
      }, this->vector(this->name()));
   }

}
}
}
