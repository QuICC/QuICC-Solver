/**
 * @file IFftOperator.cpp
 * @brief Source of the interface for a generic FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/IFftOperator.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

   void IFftOperator::cleanup()
   {
   }

   void IFftOperator::initBase() const
   {
      // Initialize FFT backend
      this->initBackend();

      // Operator specific initialization
      this->initOperator();

      // Set initialization flag
      this->mIsInitialized = true;
   }

   // anelastic overload
   void IFftOperator::initBase(std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      // Initialize FFT backend
      this->initBackendAnelastic(pF);

      // Operator specific initialization
      this->initOperator();

      // Set initialization flag
      this->mIsInitialized = true;
   }

   void IFftOperator::initOperator() const
   {
   }

   MHDFloat IFftOperator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void IFftOperator::dealias(Eigen::Ref<MatrixZ> deAliased, const Eigen::Ref<const MatrixZ>& aliased) const
   {
      std::logic_error("method not implemented!");
   }

} // namespace Fft
} // namespace Transform
} // namespace QuICC
