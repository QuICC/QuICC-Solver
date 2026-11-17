/**
 * @file IComplexOperator.cpp
 * @brief Source of the interface for a generic Complex FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/IComplexOperator.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

   void IComplexOperator::init(SharedTransformSetup spSetup) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = std::dynamic_pointer_cast<IComplexOperator::SetupType>(spSetup);

      //
      this->initBase();
   }

   void IComplexOperator::init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      throw std::logic_error("Unused interface");
   }

}
}
}
}
}
