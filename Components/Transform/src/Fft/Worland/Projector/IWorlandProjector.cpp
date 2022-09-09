/**
 * @file IWorlandProjector.cpp
 * @brief Source of the interface for a generic FFT based Worland projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Debug includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Worland/Projector/IWorlandProjector.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   IWorlandProjector::IWorlandProjector()
   {
      this->mProfileTag += "-Projector";
   }

   IWorlandProjector::~IWorlandProjector()
   {
   }

   void IWorlandProjector::initBackend() const
   {
      int lshift = 0; // operator doesn't shift l
      int extraN = 0; // no extra modes are required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void IWorlandProjector::transformBlock(MatrixZ& rOut, const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      Profiler::RegionStart<3> (this->mProfileTag + "-pre");
      this->applyPreOperator(in, isEven, useReal);
      Profiler::RegionStop<3> (this->mProfileTag + "-pre");

      Profiler::RegionStart<3> (this->mProfileTag + "-fft");
      this->mBackend.applyFft();
      Profiler::RegionStop<3> (this->mProfileTag + "-fft");

      Profiler::RegionStart<3> (this->mProfileTag + "-post");
      this->applyPostOperator(rOut, isEven, useReal);
      Profiler::RegionStop<3> (this->mProfileTag + "-post");
   }

   void IWorlandProjector::transformBlock(Matrix& rOut, const Matrix& in, const bool isEven) const
   {
      Profiler::RegionStart<3> (this->mProfileTag + "-pre");
      this->applyPreOperator(in, isEven);
      Profiler::RegionStop<3> (this->mProfileTag + "-pre");

      Profiler::RegionStart<3> (this->mProfileTag + "-fft");
      this->mBackend.applyFft();
      Profiler::RegionStop<3> (this->mProfileTag + "-fft");

      Profiler::RegionStart<3> (this->mProfileTag + "-post");
      this->applyPostOperator(rOut, isEven);
      Profiler::RegionStop<3> (this->mProfileTag + "-post");
   }

   void IWorlandProjector::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<2> fix(this->mProfileTag);

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == rOut.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true, true);
      this->transformBlock(rOut, in, true, false);
      this->transformBlock(rOut, in, false, true);
      this->transformBlock(rOut, in, false, false);
   }

   void IWorlandProjector::transform(Matrix& rOut, const Matrix& in) const
   {
      Profiler::RegionFixture<2> fix(this->mProfileTag);

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == rOut.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true);
      this->transformBlock(rOut, in, false);
   }

   int IWorlandProjector::outRows() const
   {
      return this->mspSetup->fwdSize();
   }

   int IWorlandProjector::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IWorlandProjector::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void IWorlandProjector::transform(Matrix& rOut, const MatrixZ& in) const
   {
      IWorlandOperator::transform(rOut, in);
   }

}
}
}
}
}
