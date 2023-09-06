/**
 * @file IChebyshevIntegrator.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/IChebyshevIntegrator.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Integrator {

   IChebyshevIntegrator::IChebyshevIntegrator()
   {
   }

   IChebyshevIntegrator::~IChebyshevIntegrator()
   {
   }

   void IChebyshevIntegrator::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void IChebyshevIntegrator::transform(MatrixZ& rOut, const MatrixZ& in) const
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

   void IChebyshevIntegrator::transform(Matrix& rOut, const Matrix& in) const
   {
      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() >= this->outRows());
      assert(in.cols() <= rOut.cols());

      this->mBackend.applyFft(rOut, in);
      this->applyPostOperator(rOut);
   }

   void IChebyshevIntegrator::transform(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT integrator");
   }

   void IChebyshevIntegrator::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT integrator");
   }

   int IChebyshevIntegrator::outRows() const
   {
      return this->mspSetup->specSize();
   }

   int IChebyshevIntegrator::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IChebyshevIntegrator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
}
}
