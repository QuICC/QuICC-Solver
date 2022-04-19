/**
 * @file P.cpp
 * @brief Source of the implementation of the Chebyshev P projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

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

   void P::applyPreOperator(Matrix& rOut, const Matrix& in) const
   {
      this->mBackend.input(in, true);

      this->mBackend.output(rOut);
   }

   void P::applyPostOperator(Matrix&) const
   {
   }

   void P::applyPreOperator(const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(in, useReal, true);

      this->mBackend.io();
   }

   void P::applyPostOperator(MatrixZ& rOut, const bool useReal) const
   {
      this->mBackend.output(rOut, useReal);
   }

}
}
}
}
}
}
