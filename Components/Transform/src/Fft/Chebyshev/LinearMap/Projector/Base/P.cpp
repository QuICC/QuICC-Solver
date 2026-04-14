/**
 * @file P.cpp
 * @brief Source of the implementation of the Chebyshev P projector, with linear
 * map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/P.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

void P<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const
{
   this->mBackend.input(tmp, in);
}

void P<base_t>::applyPostOperator(Eigen::Ref<Matrix>) const {}

void P<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
}

void P<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp,
   const bool useReal) const
{
   this->mBackend.output(rOut, tmp, useReal);
}

} // namespace Projector
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
