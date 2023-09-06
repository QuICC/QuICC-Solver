/**
 * @file IMixedProjector.cpp
 * @brief Source of the interface for a generic Fourier mixed projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Projector/IMixedProjector.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Projector {

   IMixedProjector::IMixedProjector()
   {
   }

   IMixedProjector::~IMixedProjector()
   {
   }

   void IMixedProjector::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void IMixedProjector::transform(Matrix& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<2> fix("IMixedProjector::transform");

      assert(this->isInitialized());
      assert(in.cols() <= rOut.cols());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() == this->outRows());

      auto& tmp = this->mBackend.getStorage();
      this->applyPreOperator(tmp, in);
      this->mBackend.applyFft(rOut, tmp);
   }

   void IMixedProjector::transform(Matrix&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Mixed FFT projector");
   }

   void IMixedProjector::transform(MatrixZ&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Mixed FFT projector");
   }

   void IMixedProjector::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Mixed FFT projector");
   }

   int IMixedProjector::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IMixedProjector::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IMixedProjector::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
      mem += this->mBackend.requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
}
}
