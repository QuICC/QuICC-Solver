/**
 * @file Spectrum.cpp
 * @brief Source of the implementation of the Chebyshev energy reductor, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Base/Spectrum.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   void Spectrum<base_t>::initOperator() const
   {
      // Check for division by 0!
      assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);
   }

   void Spectrum<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const
   {
      this->mBackend.input(tmp, in);
   }

   // I don't think we need this
   void Spectrum<base_t>::applyPostOperator(Eigen::Ref<Matrix> rOut, const Matrix& tmp) const
   {
      assert(rOut.cols() == 1);
      this->mBackend.output(rOut, tmp);
   }

   void Spectrum<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in, const bool useReal) const
   {
      this->mBackend.input(tmp, in, useReal);
   }

}
}
}
}
}
}
