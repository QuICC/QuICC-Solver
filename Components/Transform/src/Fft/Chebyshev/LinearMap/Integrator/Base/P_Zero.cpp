/**
 * @file P_Zero.cpp
 * @brief Source of the implementation of the Chebyshev P integrator, with
 * linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Base/P_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

void P_Zero<base_t>::initOperator() const
{
   if (this->mspSetup->slowSize() > 0 && this->mspSetup->slow(0) == 0)
   {
      int size =
         this->mspSetup->specSize() + std::min(0, this->mspSetup->padSize());
      this->mBackend.setMeanOperator(
         SparseMatrix(size, size).topRows(this->mspSetup->specSize()));
   }
}

void P_Zero<base_t>::applyPostOperator(Eigen::Ref<Matrix> rOut) const
{
   this->mBackend.output(rOut);
}

void P_Zero<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
}

void P_Zero<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp,
   const bool useReal) const
{
   this->mBackend.output(rOut, tmp, useReal);
}

} // namespace Integrator
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
