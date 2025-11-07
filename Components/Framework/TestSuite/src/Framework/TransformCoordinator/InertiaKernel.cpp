/**
 * @file InertiaKernel.cpp
 * @brief Source of physical space inertia kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/InertiaKernel.hpp"

// Project includes
//
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

   InertiaKernel::InertiaKernel()
      : IPhysicalKernel()
   {
   }

   std::size_t InertiaKernel::name() const
   {
      return this->mName;
   }

   void InertiaKernel::setVelocity(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

      this->mName = name;

      this->setField(name, spField);
   }

   void InertiaKernel::init(const MHDFloat inertia)
   {
      // Set scaling constants
      this->mInertia = inertia;
   }

   void InertiaKernel::setMesh(std::shared_ptr<std::vector<Array> > spMesh)
   {
      IPhysicalKernel::setMesh(spMesh);

      this->mRadius = this->mspMesh->at(0);
   }

   void InertiaKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp, FieldComponents::Physical::Id id) const
   {
      ///
      /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
      ///
      switch(id)
      {
         case(FieldComponents::Physical::R):
            std::visit([&](auto&& v){Physical::Cross<FieldComponents::Physical::THETA,FieldComponents::Physical::PHI>::set(rNLComp, v->dom(0).curl(), v->dom(0).phys(), this->mInertia);}, this->vector(this->name()));
            break;
         case(FieldComponents::Physical::THETA):
            std::visit([&](auto&& v){Physical::Cross<FieldComponents::Physical::PHI,FieldComponents::Physical::R>::set(rNLComp, v->dom(0).curl(), v->dom(0).phys(), this->mInertia);}, this->vector(this->name()));
            break;
         case(FieldComponents::Physical::PHI):
            std::visit([&](auto&& v){Physical::Cross<FieldComponents::Physical::R,FieldComponents::Physical::THETA>::set(rNLComp, v->dom(0).curl(), v->dom(0).phys(), this->mInertia);}, this->vector(this->name()));
            break;
         default:
            assert(false);
            break;
      }
   }

}
}
}
