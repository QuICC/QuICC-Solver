/**
 * @file IComplexIntegrator.cpp
 * @brief Source of the interface for a generic FFT based complex integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Complex/Integrator/IComplexIntegrator.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Complex {

namespace Integrator {

   IComplexIntegrator::IComplexIntegrator()
   {
   }

   IComplexIntegrator::~IComplexIntegrator()
   {
   }

   void IComplexIntegrator::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void IComplexIntegrator::transform(MatrixZ& rOut, const MatrixZ& in) const
   {
      Profiler::RegionFixture<2> fix("IComplexIntegrator::transform");

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() >= this->outRows());
      assert(in.cols() <= rOut.cols());

      this->mBackend.applyFft(rOut, in);
      this->applyPostOperator(rOut);
   }

   void IComplexIntegrator::transform(Matrix&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Complex FFT integrator");
   }

   void IComplexIntegrator::transform(Matrix&, const MatrixZ&) const
   {
      throw std::logic_error("Data is not compatible with Complex FFT integrator");
   }

   void IComplexIntegrator::transform(MatrixZ&, const Matrix&) const
   {
      throw std::logic_error("Data is not compatible with Complex FFT integrator");
   }

   int IComplexIntegrator::outRows() const
   {
      return this->mspSetup->specSize();
   }

   int IComplexIntegrator::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IComplexIntegrator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void IComplexIntegrator::dealias(MatrixZ& deAliased, const MatrixZ& aliased)const
   {
      int specSize = this->mspSetup->specSize();
      assert(deAliased.rows() == specSize);

      int negN = specSize/2;
      int posN = negN + (specSize%2);

      deAliased.topRows(posN) = aliased.topRows(posN);
      deAliased.bottomRows(negN) = aliased.bottomRows(negN);
   }

}
}
}
}
}
}
