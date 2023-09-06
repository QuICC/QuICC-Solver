/**
 * @file IChebyshevProjector.cpp
 * @brief Source of the interface for a generic FFT based Chebyshev projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/IChebyshevProjector.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

   IChebyshevProjector::IChebyshevProjector()
   {
   }

   IChebyshevProjector::~IChebyshevProjector()
   {
   }

   void IChebyshevProjector::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void IChebyshevProjector::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());
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

   void IChebyshevProjector::transform(Matrix& rOut, const Matrix& in) const
   {
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());
      assert(in.cols() <= rOut.cols());

      auto& tmp = this->mBackend.getStorage();
      this->applyPreOperator(tmp, in);
      this->mBackend.applyFft(rOut, tmp);
      this->applyPostOperator(rOut);
   }

   void IChebyshevProjector::transform(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT projector");
   }

   void IChebyshevProjector::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Chebyshev FFT projector");
   }

   int IChebyshevProjector::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IChebyshevProjector::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IChebyshevProjector::requiredStorage() const
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
