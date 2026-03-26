/**
 * @file ILinearMapRadialPower.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev radial power reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/ILinearMapRadialPower.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   void ILinearMapRadialPower::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void ILinearMapRadialPower::transform(Matrix& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
      this->applyPreOperator(tmpIn, in, true);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, true);
      this->applyPreOperator(tmpIn, in, false);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, false);
      this->applyPostOperator(rOut, tmpSquare);
   }

   void ILinearMapRadialPower::transform(Matrix& rOut, const Matrix& in) const
   {
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());

      auto& tmpIn = this->mBackend.getStorage(StorageKind::in);
      auto& tmpOut = this->mBackend.getStorage(StorageKind::out);
      auto& tmpSquare = this->mBackend.getStorage(StorageKind::mid);
      this->applyPreOperator(tmpIn, in);
      this->mBackend.applyFft(tmpOut, tmpIn);
      this->mBackend.square(tmpSquare, tmpOut, true);
      this->applyPostOperator(rOut, tmpSquare);
   }

   MHDFloat ILinearMapRadialPower::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   int ILinearMapRadialPower::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int ILinearMapRadialPower::outCols() const
   {
      return this->mspSetup->blockSize();
   }

}
}
}
}
}
}
