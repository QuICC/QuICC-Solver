/**
 * @file MixedProjector.cpp
 * @brief Source of the interface for a generic API for mixed projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/MixedProjector.hpp"

// Project includes
//
#if defined QUICC_FFT_MIXED_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/MixedProjector.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_MIXED_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/MixedProjector.hpp"
   #define BACKENDIMPL CuFft
#endif
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct MixedProjector::BackendImpl: public BACKENDIMPL::MixedProjector
   {
      BackendImpl() = default;
   };

   MixedProjector::MixedProjector()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   MixedProjector::~MixedProjector()
   {
   }

   void MixedProjector::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void MixedProjector::input(MatrixZ& out, const MatrixZ& in) const
   {
      this->mpImpl->input(out, in);
   }

   void MixedProjector::inputDiff(MatrixZ& out, const MatrixZ& in, const int order, const MHDFloat scale) const
   {
      this->mpImpl->inputDiff(out, in, order, scale);
   }

   void MixedProjector::applyFft(Matrix& phys, const MatrixZ& mods) const
   {
      this->mpImpl->applyFft(phys, mods);
   }

   MHDFloat MixedProjector::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
      mem += this->mpImpl->requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   MatrixZ& MixedProjector::getStorage() const
   {
      return this->mpImpl->getStorage();
   }

}
}
}
}
