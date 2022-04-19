/**
 * @file IPhysicalKernel.cpp
 * @brief Source of building block for the implementation of a physical kernel
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

// Project includes
//

namespace QuICC {

namespace Physical {

namespace Kernel {

   IPhysicalKernel::IPhysicalKernel()
   {
   }

   IPhysicalKernel::~IPhysicalKernel()
   {
   }

   void IPhysicalKernel::setMesh(std::shared_ptr<std::vector<Array> > spMesh)
   {
      this->mspMesh = spMesh;
   }

   void IPhysicalKernel::setResolution(SharedCResolution spRes)
   {
      this->mspRes = spRes;
   }

   SharedCResolution IPhysicalKernel::spRes() const
   {
      assert(this->mspRes);

      return this->mspRes;
   }

   const Resolution& IPhysicalKernel::res() const
   {
      assert(this->mspRes);

      return *this->mspRes;
   }

   void IPhysicalKernel::setField(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 0);

      this->mScalars.insert(std::make_pair(name, spField));
   }

   void IPhysicalKernel::setField(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 0);

      this->mVectors.insert(std::make_pair(name, spField));
   }

   const Framework::Selector::VariantSharedScalarVariable& IPhysicalKernel::scalar(std::size_t name) const
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 1);

      return this->mScalars.find(name)->second;
   }

   const Framework::Selector::VariantSharedVectorVariable& IPhysicalKernel::vector(std::size_t name) const
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 1);

      return this->mVectors.find(name)->second;
   }

   void IPhysicalKernel::compute(Framework::Selector::ComplexScalarField&, FieldComponents::Physical::Id) const
   {
      throw std::logic_error("Kernel does not implement complex valued computation");
   }

}
}
}
