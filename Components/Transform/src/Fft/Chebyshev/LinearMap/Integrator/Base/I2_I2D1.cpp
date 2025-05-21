/**
 * @file I2_I2D1.cpp
 * @brief Source of the implementation of the Chebyshev I^2 of P integrator, but
 * 0 mode is I^2 of D integrator, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I2D1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/Base/I2_I2D1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

void I2_I2D1<base_t>::initOperator() const
{
   int size =
      this->mspSetup->specSize() + std::min(2, this->mspSetup->padSize());
   ::QuICC::SparseSM::Chebyshev::LinearMap::I2 op(size, size,
      this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.setSpectralOperator(
      op.mat().topRows(this->mspSetup->specSize()));

   if (this->mspSetup->slowSize() > 0 && this->mspSetup->slow(0) == 0)
   {
      size =
         this->mspSetup->specSize() + std::min(1, this->mspSetup->padSize());
      ::QuICC::SparseSM::Chebyshev::LinearMap::I2D1 meanOp(size, size,
         this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setMeanOperator(
         meanOp.mat().topRows(this->mspSetup->specSize()));
   }
}

void I2_I2D1<base_t>::applyPostOperator(Matrix& rOut) const
{
   this->mBackend.outputSpectral(rOut);
}

void I2_I2D1<base_t>::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
}

void I2_I2D1<base_t>::applyPostOperator(MatrixZ& rOut, const Matrix& tmp,
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
