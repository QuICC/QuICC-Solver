/**
 * @file IMixedOperator.cpp
 * @brief Source of the interface for a generic Mixed FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/IMixedOperator.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

   IMixedOperator::IMixedOperator()
   {
   }

   IMixedOperator::~IMixedOperator()
   {
   }

   void IMixedOperator::init(SharedTransformSetup spSetup) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = std::dynamic_pointer_cast<IMixedOperator::SetupType>(spSetup);

      //
      this->initBase();
   }

   void IMixedOperator::init(SharedTransformSetup spSetup, const Internal::Array& igrid, const Internal::Array& iweights) const
   {
      throw std::logic_error("Unused interface");
   }

}
}
}
}
}
