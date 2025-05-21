/**
 * @file Y1_Zero.cpp
 * @brief Source of the implementation of the  Chebyshev Y integrator, with
 * linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Base/Y1_Zero.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

void Y1_Zero<base_t>::initOperator() const
{
   int size =
      this->mspSetup->specSize() + std::min(1, this->mspSetup->padSize());
   ::QuICC::SparseSM::Chebyshev::LinearMap::Y1 op(size, size,
      this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.setSpectralOperator(
      op.mat().topRows(this->mspSetup->specSize()));

   if (this->mspSetup->slowSize() > 0 && this->mspSetup->slow(0) == 0)
   {
      this->mBackend.setMeanOperator(
         SparseMatrix(size, size).topRows(this->mspSetup->specSize()));
   }
}

void Y1_Zero<base_t>::applyPostOperator(Matrix& rOut) const
{
   this->mBackend.outputSpectral(rOut);
}

void Y1_Zero<base_t>::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
}

void Y1_Zero<base_t>::applyPostOperator(MatrixZ& rOut, const Matrix& tmp,
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
