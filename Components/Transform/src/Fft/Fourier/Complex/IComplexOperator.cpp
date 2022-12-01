/**
 * @file IComplexOperator.cpp
 * @brief Source of the interface for a generic Complex FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/IComplexOperator.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

   IComplexOperator::IComplexOperator()
   {
   }

   IComplexOperator::~IComplexOperator()
   {
   }

   void IComplexOperator::init(SharedTransformSetup spSetup) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = std::dynamic_pointer_cast<IComplexOperator::SetupType>(spSetup);

      //
      this->initBase();
   }

   void IComplexOperator::init(SharedTransformSetup spSetup, const internal::Array& igrid, const internal::Array& iweights) const
   {
      throw std::logic_error("Unused interface");
   }


}
}
}
}
}
