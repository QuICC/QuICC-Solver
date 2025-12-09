/**
 * @file IChebyshevOperator.cpp
 * @brief Source of the interface for a generic Chebyshev FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

   void IChebyshevOperator::cleanup()
   {
   }

   void IChebyshevOperator::init(SharedTransformSetup spSetup) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = std::dynamic_pointer_cast<IChebyshevOperator::SetupType>(spSetup);

      //
      this->initBase();
   }

   // anelastic overload
   void IChebyshevOperator::init(SharedTransformSetup spSetup, std::shared_ptr<QuICC::DenseSM::Chebyshev::LinearMap::RadialTorPolFunction> pF) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = std::dynamic_pointer_cast<IChebyshevOperator::SetupType>(spSetup);

      //
      this->initBase(pF);
   }

   void IChebyshevOperator::init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      throw std::logic_error("Unused interface");
   }

   void IChebyshevOperator::initBase() const
   {
      // Initialize FFT backend
      this->initBackend();

      // Operator specific initialization
      this->initOperator();

      // Set initialization flag
      this->mIsInitialized = true;
   }

   void IChebyshevOperator::initOperator() const
   {
   }

   MHDFloat IChebyshevOperator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
