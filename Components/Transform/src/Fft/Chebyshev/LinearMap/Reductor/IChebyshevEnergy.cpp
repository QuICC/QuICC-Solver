/**
 * @file IChebyshevEnergy.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev energy reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/IChebyshevEnergy.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

   IChebyshevEnergy::IChebyshevEnergy()
   {
   }

   IChebyshevEnergy::~IChebyshevEnergy()
   {
   }

   void IChebyshevEnergy::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void IChebyshevEnergy::transform(Matrix& rOut, const MatrixZ& in) const
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
      this->mBackend.applyFwdFft(tmpOut, tmpSquare);
      this->applyPostOperator(rOut, tmpOut);
   }

   void IChebyshevEnergy::transform(Matrix& rOut, const Matrix& in) const
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
      this->mBackend.applyFwdFft(tmpOut, tmpSquare);
      this->applyPostOperator(rOut, tmpOut);
   }

   void IChebyshevEnergy::transform(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT energy reductor");
   }

   void IChebyshevEnergy::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT energy reductor");
   }

   MHDFloat IChebyshevEnergy::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   int IChebyshevEnergy::outRows() const
   {
      return this->mspSetup->blockSize();
   }

   int IChebyshevEnergy::outCols() const
   {
      return 1;
   }

}
}
}
}
}
}
