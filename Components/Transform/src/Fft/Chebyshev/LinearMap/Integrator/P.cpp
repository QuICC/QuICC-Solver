/**
 * @file P.cpp
 * @brief Source of the implementation of the Chebyshev P integrator, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   P::P()
   {
   }

   P::~P()
   {
   }

   void P::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.output(rOut);
   }

   void P::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

   void P::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.output(rOut, tmp, useReal);
   }

}
}
}
}
}
}
