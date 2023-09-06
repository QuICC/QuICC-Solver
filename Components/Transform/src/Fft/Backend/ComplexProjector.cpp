/**
 * @file ComplexProjector.cpp
 * @brief Source of the interface for a generic API for complex projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/ComplexProjector.hpp"

// Project includes
//
#if defined QUICC_FFT_COMPLEX_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/ComplexProjector.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_COMPLEX_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/ComplexProjector.hpp"
   #define BACKENDIMPL CuFft
#endif
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct ComplexProjector::BackendImpl: public BACKENDIMPL::ComplexProjector
   {
      BackendImpl() = default;
   };

   ComplexProjector::ComplexProjector()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   ComplexProjector::~ComplexProjector()
   {
   }

   void ComplexProjector::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void ComplexProjector::initMeanBlocks(const MatrixI& idBlocks) const
   {
      this->mpImpl->initMeanBlocks(idBlocks);
   }

   void ComplexProjector::applyFft(MatrixZ& phys, const MatrixZ& mods) const
   {
      this->mpImpl->applyFft(phys, mods);
   }

   void ComplexProjector::input(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mpImpl->input(tmp, in);
   }

   void ComplexProjector::inputMean(MatrixZ& tmp, const MatrixZ& in) const
   {
      this->mpImpl->inputMean(tmp, in);
   }

   void ComplexProjector::inputDiff(MatrixZ& tmp, const MatrixZ& in, const int order, const MHDFloat scale) const
   {
      this->mpImpl->inputDiff(tmp, in, order, scale);
   }

   void ComplexProjector::inputDiff2D(MatrixZ& tmp, const MatrixZ& in, const std::vector<std::pair<int,int> >& orders, const MHDFloat scale, const MatrixI& idBlocks) const
   {
      this->mpImpl->inputDiff2D(tmp, in, orders, scale, idBlocks);
   }

   MHDFloat ComplexProjector::requiredStorage() const
   {
      MHDFloat mem = 0.0;

#ifdef QUICC_STORAGEPROFILE
      mem += static_cast<MHDFloat>(Debug::MemorySize<int>::BYTES);
      mem += this->mpImpl->requiredStorage();
#endif // QUICC_STORAGEPROFILE

      return mem;
   }

   MatrixZ& ComplexProjector::getStorage() const
   {
      return this->mpImpl->getStorage();
   }

}
}
}
}
