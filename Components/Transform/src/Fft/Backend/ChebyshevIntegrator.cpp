/**
 * @file ChebyshevIntegrator.cpp
 * @brief Source of the interface for generic aPI for Chebyshev integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/ChebyshevIntegrator.hpp"

// Project includes
//
#if defined QUICC_FFT_CHEBYSHEV_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/ChebyshevIntegrator.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_CHEBYSHEV_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/ChebyshevIntegrator.hpp"
   #define BACKENDIMPL CuFft
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct ChebyshevIntegrator::BackendImpl: public BACKENDIMPL::ChebyshevIntegrator
   {
      BackendImpl() = default;
   };

   ChebyshevIntegrator::ChebyshevIntegrator()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   ChebyshevIntegrator::~ChebyshevIntegrator()
   {
   }

   void ChebyshevIntegrator::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void ChebyshevIntegrator::input(Matrix& tmp, const MatrixZ& in, const bool useReal) const
   {
      this->mpImpl->input(tmp, in, useReal);
   }

   void ChebyshevIntegrator::applyFft(Matrix& mods, const Matrix& phys) const
   {
      this->mpImpl->applyFft(mods, phys);
   }

   void ChebyshevIntegrator::setSpectralOperator(const SparseMatrix& mat) const
   {
      this->mpImpl->setSpectralOperator(mat);
   }

   void ChebyshevIntegrator::setMeanOperator(const SparseMatrix& mat) const
   {
      this->mpImpl->setMeanOperator(mat);
   }

   void ChebyshevIntegrator::output(Matrix& rOut) const
   {
      this->mpImpl->output(rOut);
   }

   void ChebyshevIntegrator::outputSpectral(Matrix& rOut) const
   {
      this->mpImpl->outputSpectral(rOut);
   }

   void ChebyshevIntegrator::output(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mpImpl->output(rOut, tmp, useReal);
   }

   void ChebyshevIntegrator::outputSpectral(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mpImpl->outputSpectral(rOut, tmp, useReal);
   }

   Matrix& ChebyshevIntegrator::getStorage(const StorageKind kind) const
   {
      return this->mpImpl->getStorage(kind);
   }

}
}
}
}
