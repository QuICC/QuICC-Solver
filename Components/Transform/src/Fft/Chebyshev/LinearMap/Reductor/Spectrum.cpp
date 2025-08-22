/**
 * @file Spectrum.cpp
 * @brief Source of the implementation of the Chebyshev energy reductor, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Spectrum.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   Spectrum::Spectrum()
   {
   }

   Spectrum::~Spectrum()
   {
   }

   void Spectrum::applyPreOperator(Matrix& tmp, const Matrix& in) const
   {
      this->mBackend.input(tmp, in);
   }

   void Spectrum::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
   {
      assert(rOut.cols() == 1);
      this->mBackend.output(rOut, tmp);
   }

   void Spectrum::applyPreOperator(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

}
}
}
}
}
}
