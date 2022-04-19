/**
 * @file IChebyshevOperator.cpp
 * @brief Source of the interface for a generic Chebyshev FFT based operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/IChebyshevOperator.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

   IChebyshevOperator::IChebyshevOperator()
   {
   }

   IChebyshevOperator::~IChebyshevOperator()
   {
   }

   void IChebyshevOperator::init(IChebyshevOperator::SharedSetupType spSetup) const
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      //
      this->initBase();
   }

}
}
}
}
