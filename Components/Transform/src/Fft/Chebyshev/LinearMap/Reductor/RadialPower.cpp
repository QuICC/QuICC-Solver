/**
 * @file RadialPower.cpp
 * @brief Source of the implementation of the Chebyshev radial power reductor,
 * with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/RadialPower.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

void RadialPower::applyPreOperator(Matrix& tmp, const Matrix& in) const
{
   this->mBackend.input(tmp, in);
}

void RadialPower::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
{
   assert(rOut.cols() == 1);
   this->mBackend.outputGrid(rOut, tmp);
}

void RadialPower::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
