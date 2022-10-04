/**
 * @file Energy.cpp
 * @brief Source of the implementation of the Chebyshev energy reductor, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Energy.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   Energy::Energy()
   {
   }

   Energy::~Energy()
   {
   }

   void Energy::applyPreOperator(Matrix& tmp, const Matrix& in) const
   {
      this->mBackend.input(tmp, in);
   }

   void Energy::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
   {
      assert(rOut.cols() == 1);
      this->mBackend.output(rOut, tmp);
   }

   void Energy::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

}
}
}
}
}
}
