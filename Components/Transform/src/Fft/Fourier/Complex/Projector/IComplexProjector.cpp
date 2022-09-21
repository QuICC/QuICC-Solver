/**
 * @file IComplexProjector.cpp
 * @brief Source of the interface for a generic Fourier complex projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Complex/Projector/IComplexProjector.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Projector {

   IComplexProjector::IComplexProjector()
   {
   }

   IComplexProjector::~IComplexProjector()
   {
   }

   void IComplexProjector::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void IComplexProjector::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      assert(this->isInitialized());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());
      assert(in.cols() <= rOut.cols());

      auto& tmp = this->mBackend.getStorage();
      this->applyPreOperator(tmp, in);
      this->mBackend.applyFft(rOut, tmp);
   }

   void IComplexProjector::transform(Matrix&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Complex FFT projector");
   }

   void IComplexProjector::transform(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Complex FFT projector");
   }

   void IComplexProjector::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Complex FFT projector");
   }

   int IComplexProjector::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IComplexProjector::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IComplexProjector::requiredStorage() const
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
