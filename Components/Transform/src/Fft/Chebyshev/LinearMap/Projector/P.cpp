/**
 * @file P.cpp
 * @brief Source of the implementation of the Chebyshev P projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   P::P()
   {
   }

   P::~P()
   {
   }

   void P::applyPreOperator(Matrix& tmp, const Matrix& in) const
   {
      this->mBackend.input(tmp, in);
   }

   void P::applyPostOperator(Matrix&) const
   {
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
