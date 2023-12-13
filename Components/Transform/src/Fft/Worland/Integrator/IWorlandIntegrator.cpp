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
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Worland {

namespace Integrator {

   IWorlandIntegrator::IWorlandIntegrator()
   {
      this->mProfileTag += "-Integrator";
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
      Profiler::RegionStart<5> (this->mProfileTag + "-pre");
      this->applyPreOperator(in, isEven, useReal);
      Profiler::RegionStop<5> (this->mProfileTag + "-pre");

      Profiler::RegionStart<5> (this->mProfileTag + "-fft");
      this->mBackend.applyFft();
      Profiler::RegionStop<5> (this->mProfileTag + "-fft");

      Profiler::RegionStart<5> (this->mProfileTag + "-post");
      this->applyPostOperator(rOut, isEven, useReal);
      Profiler::RegionStop<5> (this->mProfileTag + "-post");
   }

   void IWorlandIntegrator::transformBlock(Matrix& rOut, const Matrix& in, const bool isEven) const
   {
      Profiler::RegionStart<5> (this->mProfileTag + "-pre");
      this->applyPreOperator(in, isEven);
      Profiler::RegionStop<5> (this->mProfileTag + "-pre");

      Profiler::RegionStart<5> (this->mProfileTag + "-fft");
      this->mBackend.applyFft();
      Profiler::RegionStop<5> (this->mProfileTag + "-fft");

      Profiler::RegionStart<5> (this->mProfileTag + "-post");
      this->applyPostOperator(rOut, isEven);
      Profiler::RegionStop<5> (this->mProfileTag + "-post");
   }

   void IWorlandIntegrator::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<2> fix(this->mProfileTag + "::transform");

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true, true);
      this->transformBlock(rOut, in, true, false);
      this->transformBlock(rOut, in, false, true);
      this->transformBlock(rOut, in, false, false);
   }

   void IWorlandIntegrator::transform(Matrix& rOut, const Matrix& in) const
   {
      Profiler::RegionFixture<2> fix(this->mProfileTag + "::transform");

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true);
      this->transformBlock(rOut, in, false);
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

   void IWorlandIntegrator::transform(MatrixZ& rOut, const Matrix& in) const
   {
      IWorlandOperator::transform(rOut, in);
   }

}
}
}
}
}
