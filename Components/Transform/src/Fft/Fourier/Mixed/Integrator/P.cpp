/**
 * @file P.cpp
 * @brief Source of the implementation of the Fourier mixed P integrator
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   P::P()
   {
   }

   P::~P()
   {
   }

   void P::applyPreOperator(MatrixZ& rOut, const Matrix& in) const
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
