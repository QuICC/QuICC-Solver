/**
 * @file P.cpp
 * @brief Source of the implementation of the Chebyshev P integrator, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Base/P.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   void P<base_t>::applyPostOperator(Matrix& rOut) const
   {
      this->mBackend.output(rOut);
   }

   void P<base_t>::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

   void P<base_t>::applyPostOperator(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mBackend.output(rOut, tmp, useReal);
   }

}
}
}
}
}
}
