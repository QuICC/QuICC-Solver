/**
 * @file IALegendreIntegrator.cpp
 * @brief Source of the interface for a generic FFT based ALegendre integrator
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
#include "QuICC/Transform/Fft/ALegendre/Integrator/IALegendreIntegrator.hpp"

// Project includes
//
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace ALegendre {

namespace Integrator {

   IALegendreIntegrator::IALegendreIntegrator()
   {
      this->mProfileTag += "-Integrator";
   }

   IALegendreIntegrator::~IALegendreIntegrator()
   {
   }

   void IALegendreIntegrator::initBackend() const
   {
      int lshift = 0; // operator doesn't shift l
      int extraN = 0; // no extra modes are required
      this->mBackend.init(*this->mspSetup, lshift, extraN);
   }

   void IALegendreIntegrator::transformBlock(MatrixZ& rOut, const MatrixZ& in, const bool isEven, const bool useReal) const
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

   void IALegendreIntegrator::transformBlock(Matrix& rOut, const Matrix& in, const bool isEven) const
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

   void IALegendreIntegrator::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<2> fix(this->mProfileTag + "::transform");

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(in.cols() <= rOut.cols());
      /* for (int i = 0; i < 5; i++)
      {
         for (int j = 0; j < 5; j++)
         {
            printf("%.2e %.2e | ", in(j, i).real(),
               in(j, i).imag());
         }
         printf("aa\n");
      }*/
#ifdef QUICC_USE_PFSOLVE
      this->transformBlock(rOut, in, true, true);
#else
      this->transformBlock(rOut, in, true, true);
      this->transformBlock(rOut, in, true, false);
      this->transformBlock(rOut, in, false, true);
      this->transformBlock(rOut, in, false, false);
#endif
      /* for (int i = 0; i < 5; i++)
      {
         for (int j = 0; j < 5; j++)
         {
            printf("%.2e %.2e | ", rOut(j, i).real(),
               rOut(j, i).imag());
         }
         printf("bb\n");
      }*/
   }

   void IALegendreIntegrator::transform(Matrix& rOut, const Matrix& in) const
   {
      Profiler::RegionFixture<2> fix(this->mProfileTag + "::transform");

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(in.cols() <= rOut.cols());

      this->transformBlock(rOut, in, true);
      this->transformBlock(rOut, in, false);
   }

   int IALegendreIntegrator::outRows() const
   {
      return this->mspSetup->fastSize(0);
   }

   int IALegendreIntegrator::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IALegendreIntegrator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void IALegendreIntegrator::transform(Matrix& rOut, const MatrixZ& in) const
   {
      IALegendreOperator::transform(rOut, in);
   }

   void IALegendreIntegrator::transform(MatrixZ& rOut, const Matrix& in) const
   {
      IALegendreOperator::transform(rOut, in);
   }

}
}
}
}
}
