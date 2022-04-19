/**
 * @file MixedIntegrator.cpp
 * @brief Source of the interface for a generic API for mixed integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/MixedIntegrator.hpp"

// Project includes
//
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

   MixedIntegrator::~MixedIntegrator()
   {
   }

   void MixedIntegrator::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void MixedIntegrator::output(MatrixZ& rOut) const
   {
      this->mpImpl->output(rOut);
   }

   void MixedIntegrator::outputDiff(MatrixZ& rOut, const int order, const MHDFloat scale) const
   {
      this->mpImpl->outputDiff(rOut, order, scale);
   }

   void MixedIntegrator::outputDiff(MatrixZ& rOut, const int order, const MHDFloat scale, const std::map<int,MHDComplex>& mod) const
   {
      this->mpImpl->outputDiff(rOut, order, scale, mod);
   }

   void MixedIntegrator::io(MatrixZ& rOut, const Matrix& in) const
   {
      this->mpImpl->io(rOut, in);
   }

   void MixedIntegrator::io(MHDComplex* out, const MHDFloat* in) const
   {
      this->mpImpl->io(out, in);
   }

   void MixedIntegrator::applyFft() const
   {
      this->mpImpl->applyFft();
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
