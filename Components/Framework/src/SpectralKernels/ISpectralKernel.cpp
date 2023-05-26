/**
 * @file ISpectralKernel.cpp
 * @brief Source of building block for the implementation of a spectral kernel
 */

// System includes
//

// Project includes
//
#include "QuICC/SpectralKernels/ISpectralKernel.hpp"

namespace QuICC {

namespace Spectral {

namespace Kernel {

   ISpectralKernel::ISpectralKernel(const bool isComplex)
      : mIsComplex(isComplex)
   {
   }

   void ISpectralKernel::setResolution(SharedCResolution spRes)
   {
      this->mspRes = spRes;
   }

   SharedCResolution ISpectralKernel::spRes() const
   {
      assert(this->mspRes);

      return this->mspRes;
   }

   const Resolution& ISpectralKernel::res() const
   {
      assert(this->mspRes);

      return *this->mspRes;
   }

   void ISpectralKernel::apply(const std::size_t)
   {
   }

   void ISpectralKernel::setField(std::size_t name, Framework::Selector::VariantSharedScalarVariable spField)
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 0);

      this->mScalars.insert(std::make_pair(name, spField));
   }

   void ISpectralKernel::setField(std::size_t name, Framework::Selector::VariantSharedVectorVariable spField)
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 0);

      this->mVectors.insert(std::make_pair(name, spField));
   }

   const Framework::Selector::VariantSharedScalarVariable& ISpectralKernel::scalar(std::size_t name) const
   {
      // Safety assertion
      assert(this->mScalars.count(name) == 1);

      return this->mScalars.find(name)->second;
   }

   const Framework::Selector::VariantSharedVectorVariable& ISpectralKernel::vector(std::size_t name) const
   {
      // Safety assertion
      assert(this->mVectors.count(name) == 1);

      return this->mVectors.find(name)->second;
   }

} // Kernel
} // Spectral
} // QuICC
