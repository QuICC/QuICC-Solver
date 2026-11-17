/**
 * @file IMixedIntegrator.cpp
 * @brief Source of the interface for a generic FFT based mixed integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Fourier/Mixed/Integrator/IMixedIntegrator.hpp"
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Fourier {

namespace Mixed {

namespace Integrator {

   void IMixedIntegrator::initBackend() const
   {
      this->mBackend.init(*this->mspSetup);
   }

   void IMixedIntegrator::transform(Eigen::Ref<MatrixZ> rOut, const Eigen::Ref<const Matrix>& in) const
   {
      Profiler::RegionFixture<2> fix("IMixedIntegrator::transform");

      assert(this->isInitialized());
      assert(this->mspSetup->fwdSize() == in.rows());
      assert(rOut.cols() == this->outCols());
      assert(rOut.rows() >= this->outRows());
      assert(in.cols() <= rOut.cols());

      this->mBackend.applyFft(rOut, in);
      this->applyPostOperator(rOut);
   }

   void IMixedIntegrator::transform(Eigen::Ref<Matrix>, const Eigen::Ref<const Matrix>&) const
   {
      throw std::logic_error("Data is not compatible with Mixed FFT integrator");
   }

   void IMixedIntegrator::transform(Eigen::Ref<Matrix>, const Eigen::Ref<const MatrixZ>&) const
   {
      throw std::logic_error("Data is not compatible with Mixed FFT integrator");
   }

   void IMixedIntegrator::transform(Eigen::Ref<MatrixZ>, const Eigen::Ref<const MatrixZ>&) const
   {
      throw std::logic_error("Data is not compatible with Mixed FFT integrator");
   }

   int IMixedIntegrator::outRows() const
   {
      return this->mspSetup->specSize();
   }

   int IMixedIntegrator::outCols() const
   {
      return this->mspSetup->blockSize();
   }

   MHDFloat IMixedIntegrator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   void IMixedIntegrator::dealias(Eigen::Ref<MatrixZ> deAliased, const Eigen::Ref<const MatrixZ>& aliased)const
   {
      int specSize = this->mspSetup->specSize();
      assert(deAliased.rows() == specSize);

      deAliased = aliased.topRows(specSize);
   }


}
}
}
}
}
}
