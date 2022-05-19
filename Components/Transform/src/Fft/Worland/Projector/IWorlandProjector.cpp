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
#include "QuICC/Debug/Profiler/BreakPoint.hpp"
#include "QuICC/Transform/Fft/Worland/Projector/IWorlandProjector.hpp"

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Projector {

   IWorlandProjector::IWorlandProjector()
   {
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
      ProfilerMacro_start_lvl3(this->mProfileId, 1);
      this->applyPreOperator(in, isEven, useReal);
      ProfilerMacro_stop_lvl3(this->mProfileId, 1);

      ProfilerMacro_start_lvl3(this->mProfileId, 2);
      this->mBackend.applyFft();
      ProfilerMacro_stop_lvl3(this->mProfileId, 2);

      ProfilerMacro_start_lvl3(this->mProfileId, 3);
      this->applyPostOperator(rOut, isEven, useReal);
      ProfilerMacro_stop_lvl3(this->mProfileId, 3);
   }

   void IWorlandProjector::transformBlock(Matrix& rOut, const Matrix& in, const bool isEven) const
   {
      ProfilerMacro_start_lvl3(this->mProfileId, 1);
      this->applyPreOperator(in, isEven);
      ProfilerMacro_stop_lvl3(this->mProfileId, 1);

      ProfilerMacro_start_lvl3(this->mProfileId, 2);
      this->mBackend.applyFft();
      ProfilerMacro_stop_lvl3(this->mProfileId, 2);

      ProfilerMacro_start_lvl3(this->mProfileId, 3);
      this->applyPostOperator(rOut, isEven);
      ProfilerMacro_stop_lvl3(this->mProfileId, 3);
   }

   void IWorlandProjector::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      ProfilerMacro_start(Debug::Profiler::WORLANDTRA);
      ProfilerMacro_start(Debug::Profiler::WORLANDPROJ);
      ProfilerMacro_start(this->mProfileId);

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == rOut.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true, true);
      this->transformBlock(rOut, in, true, false);
      this->transformBlock(rOut, in, false, true);
      this->transformBlock(rOut, in, false, false);

      ProfilerMacro_stop(this->mProfileId);
      ProfilerMacro_stop(Debug::Profiler::WORLANDPROJ);
      ProfilerMacro_stop(Debug::Profiler::WORLANDTRA);
   }

   void IWorlandProjector::transform(Matrix& rOut, const Matrix& in) const
   {
      ProfilerMacro_start(Debug::Profiler::WORLANDTRA);
      ProfilerMacro_start(Debug::Profiler::WORLANDPROJ);
      ProfilerMacro_start(this->mProfileId);

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == rOut.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true);
      this->transformBlock(rOut, in, false);

      ProfilerMacro_stop(this->mProfileId);
      ProfilerMacro_stop(Debug::Profiler::WORLANDPROJ);
      ProfilerMacro_stop(Debug::Profiler::WORLANDTRA);
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
