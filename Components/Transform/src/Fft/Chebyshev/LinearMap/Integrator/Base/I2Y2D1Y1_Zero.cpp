/**
 * @file I2Y2D1Y1_Zero.cpp
 * @brief Source of the implementation of the  Chebyshev I2R3 of 1/R D R
 * integrator, but 0 mode is zeroed, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2Y2D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Base/I2Y2D1Y1_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

void I2Y2D1Y1_Zero<base_t>::initOperator() const
{
   int size =
      this->mspSetup->specSize() + std::min(1, this->mspSetup->padSize());
   ::QuICC::SparseSM::Chebyshev::LinearMap::I2Y2D1Y1 op(size, size,
      this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.setSpectralOperator(
      op.mat().topRows(this->mspSetup->specSize()));

   if (this->mspSetup->slowSize() > 0 && this->mspSetup->slow(0) == 0)
   {
      this->mBackend.setMeanOperator(
         SparseMatrix(size, size).topRows(this->mspSetup->specSize()));
   }
}

void I2Y2D1Y1_Zero<base_t>::applyPostOperator(Eigen::Ref<Matrix> rOut) const
{
   this->mBackend.outputSpectral(rOut);
}

void I2Y2D1Y1_Zero<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
}

void I2Y2D1Y1_Zero<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp,
   const bool useReal) const
{
   this->mBackend.outputSpectral(rOut, tmp, useReal);
}

} // namespace Integrator
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
