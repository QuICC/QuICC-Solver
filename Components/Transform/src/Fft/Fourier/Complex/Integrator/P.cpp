/**
 * @file P.cpp
 * @brief Source of the implementation of the Fourier complex P integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   P::P()
   {
   }

   P::~P()
   {
   }

   void P::applyPreOperator(MatrixZ& rOut, const MatrixZ& in) const
   {
      this->mBackend.io(rOut, in);
   }

   void P::applyPostOperator(MatrixZ& rOut) const
   {
      this->mBackend.output(rOut);
   }

}
}
}
}
}
}
