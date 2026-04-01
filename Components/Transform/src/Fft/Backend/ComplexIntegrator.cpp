/**
 * @file ComplexIntegrator.cpp
 * @brief Source of the interface for a generic API for complex integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Backend/ComplexIntegrator.hpp"
#if defined QUICC_FFT_COMPLEX_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/ComplexIntegrator.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_COMPLEX_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/ComplexIntegrator.hpp"
   #define BACKENDIMPL CuFft
#endif
#include "QuICC/Debug/StorageProfiler/MemorySize.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct ComplexIntegrator::BackendImpl: public BACKENDIMPL::ComplexIntegrator
   {
      BackendImpl() = default;
   };

   ComplexIntegrator::ComplexIntegrator()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   void ComplexIntegrator::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void ComplexIntegrator::initMeanBlocks(const MatrixI& idBlocks) const
   {
      this->mpImpl->initMeanBlocks(idBlocks);
   }

   void ComplexIntegrator::output(Eigen::Ref<MatrixZ> rOut) const
   {
      this->mpImpl->output(rOut);
   }

   void ComplexIntegrator::outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale) const
   {
      this->mpImpl->outputDiff(rOut, order, scale);
   }

   void ComplexIntegrator::zeroMean(Eigen::Ref<MatrixZ> rOut) const
   {
      this->mpImpl->zeroMean(rOut);
   }

   void ComplexIntegrator::outputMean(Eigen::Ref<MatrixZ> rOut) const
   {
      this->mpImpl->outputMean(rOut);
   }

   void ComplexIntegrator::applyFft(Eigen::Ref<MatrixZ> mods, const Eigen::Ref<const MatrixZ>& phys) const
   {
      this->mpImpl->applyFft(mods, phys);
   }

   void ComplexIntegrator::extractMean(const MatrixZ& rOut) const
   {
      this->mpImpl->extractMean(rOut);
   }

   void ComplexIntegrator::setMean(Eigen::Ref<MatrixZ> rOut, const MHDFloat scale) const
   {
      this->mpImpl->setMean(rOut, scale);
   }

   int ComplexIntegrator::computeDiff2D(const std::vector<std::pair<int,int> >& orders, const MHDFloat scale, const MatrixI& idBlocks, const bool isInverse) const
   {
      return this->mpImpl->computeDiff2D(orders, scale, idBlocks, isInverse);
   }

   int ComplexIntegrator::multDiff2D(const int idA, const int idB) const
   {
      return this->mpImpl->multDiff2D(idA, idB);
   }

   void ComplexIntegrator::applyDiff2D(Eigen::Ref<MatrixZ> rOut, const int id) const
   {
      this->mpImpl->applyDiff2D(rOut, id);
   }

   void ComplexIntegrator::destroyDiff2D(const int id) const
   {
      this->mpImpl->destroyDiff2D(id);
   }

   MHDFloat ComplexIntegrator::requiredStorage() const
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
