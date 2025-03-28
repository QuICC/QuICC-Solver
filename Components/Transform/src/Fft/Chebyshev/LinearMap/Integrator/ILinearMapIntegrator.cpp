/**
 * @file ILinearMapIntegrator.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/ILinearMapIntegrator.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

void ILinearMapIntegrator::initBackend() const
{
   this->mBackend.init(*this->mspSetup);
}

void ILinearMapIntegrator::transform(MatrixZ& rOut, const MatrixZ& in) const
{
   assert(this->isInitialized());
   assert(this->mspSetup->fwdSize() == in.rows());
   assert(rOut.cols() == this->outCols());
   assert(rOut.rows() >= this->outRows());
   assert(in.cols() <= rOut.cols());

   auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
   auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
   this->applyPreOperator(tmpIn, in, true);
   this->mBackend.applyFft(tmpOut, tmpIn);
   this->applyPostOperator(rOut, tmpOut, true);

   this->applyPreOperator(tmpIn, in, false);
   this->mBackend.applyFft(tmpOut, tmpIn);
   this->applyPostOperator(rOut, tmpOut, false);
}

void ILinearMapIntegrator::transform(Matrix& rOut, const Matrix& in) const
{
   assert(this->isInitialized());
   assert(this->mspSetup->fwdSize() == in.rows());
   assert(rOut.cols() == this->outCols());
   assert(rOut.rows() >= this->outRows());
   assert(in.cols() <= rOut.cols());

   this->mBackend.applyFft(rOut, in);
   this->applyPostOperator(rOut);
}

void ILinearMapIntegrator::transform(Matrix&, const MatrixZ&) const
{
   throw std::logic_error(
      "Data is not compatible with Chebyshev FFT integrator");
}

int ILinearMapIntegrator::outRows() const
{
   return this->mspSetup->specSize();
}

int ILinearMapIntegrator::outCols() const
{
   return this->mspSetup->blockSize();
}

MHDFloat ILinearMapIntegrator::requiredStorage() const
{
   MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
   mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

   return mem;
}

} // namespace Integrator
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
