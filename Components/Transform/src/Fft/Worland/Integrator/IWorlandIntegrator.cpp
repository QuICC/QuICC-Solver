/**
 * @file IWorlandIntegrator.cpp
 * @brief Source of the interface for a generic FFT based Worland integrator
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
#include "QuICC/Transform/Fft/Worland/Integrator/IWorlandIntegrator.hpp"

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   IWorlandIntegrator::IWorlandIntegrator()
   {
   }

   IWorlandIntegrator::~IWorlandIntegrator()
   {
   }

   void IWorlandIntegrator::initBackend() const
   {
      int lshift = 0; // operator doesn't shift l
      int extraN = 0; // no extra modes are required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void IWorlandIntegrator::transformBlock(MatrixZ& rOut, const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->applyPreOperator(in, isEven, useReal);

      this->mBackend.applyFft();

      this->applyPostOperator(rOut, isEven, useReal);
   }

   void IWorlandIntegrator::transformBlock(Matrix& rOut, const Matrix& in, const bool isEven) const
   {
      this->applyPreOperator(in, isEven);

      this->mBackend.applyFft();

      this->applyPostOperator(rOut, isEven);
   }

   void IWorlandIntegrator::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      ProfilerMacro_start(Debug::Profiler::WORLANDTRA);
      ProfilerMacro_start(Debug::Profiler::WORLANDINTG);
      ProfilerMacro_start(this->mProfileId);

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true, true);
      this->transformBlock(rOut, in, true, false);
      this->transformBlock(rOut, in, false, true);
      this->transformBlock(rOut, in, false, false);

      ProfilerMacro_stop(this->mProfileId);
      ProfilerMacro_stop(Debug::Profiler::WORLANDINTG);
      ProfilerMacro_stop(Debug::Profiler::WORLANDTRA);
   }

   void IWorlandIntegrator::transform(Matrix& rOut, const Matrix& in) const
   {
      ProfilerMacro_start(Debug::Profiler::WORLANDTRA);
      ProfilerMacro_start(Debug::Profiler::WORLANDINTG);
      ProfilerMacro_start(this->mProfileId);

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true);
      this->transformBlock(rOut, in, false);

      ProfilerMacro_stop(this->mProfileId);
      ProfilerMacro_stop(Debug::Profiler::WORLANDINTG);
      ProfilerMacro_stop(Debug::Profiler::WORLANDTRA);
   }

   int IWorlandIntegrator::outRows() const
   {
      return this->mspSetup->fastSize(0);
   }

   int IWorlandIntegrator::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IWorlandIntegrator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void IWorlandIntegrator::transform(Matrix& rOut, const MatrixZ& in) const
   {
      IWorlandOperator::transform(rOut, in);
   }

}
}
}
}
}
