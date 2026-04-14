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

}
}
}
}
