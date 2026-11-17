/**
 * @file MixedIntegrator.cpp
 * @brief Source of the interface for a generic API for mixed integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Backend/MixedIntegrator.hpp"
#if defined QUICC_FFT_MIXED_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/MixedIntegrator.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_MIXED_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/MixedIntegrator.hpp"
   #define BACKENDIMPL CuFft
#endif
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct MixedIntegrator::BackendImpl: public BACKENDIMPL::MixedIntegrator
   {
      BackendImpl() = default;
   };

   MixedIntegrator::MixedIntegrator()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   void MixedIntegrator::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void MixedIntegrator::output(Eigen::Ref<MatrixZ> rOut) const
   {
      this->mpImpl->output(rOut);
   }

   void MixedIntegrator::outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale) const
   {
      this->mpImpl->outputDiff(rOut, order, scale);
   }

   void MixedIntegrator::outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale, const std::map<int,MHDComplex>& mod) const
   {
      this->mpImpl->outputDiff(rOut, order, scale, mod);
   }

   void MixedIntegrator::applyFft(Eigen::Ref<MatrixZ> mods, const Eigen::Ref<const Matrix>& phys) const
   {
      this->mpImpl->applyFft(mods, phys);
   }

   MHDFloat MixedIntegrator::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
      mem += this->mpImpl->requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

}
}
}
}
