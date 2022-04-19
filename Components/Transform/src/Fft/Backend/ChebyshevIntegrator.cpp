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

   void ChebyshevIntegrator::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpImpl->io(out, in);
   }

   void ChebyshevIntegrator::io(Matrix& rOut, const Matrix& in) const
   {
      this->mpImpl->io(rOut, in);
   }

   void ChebyshevIntegrator::input(const MatrixZ& in, const bool useReal) const
   {
      this->mpImpl->input(in, useReal);
   }

   void ChebyshevIntegrator::applyFft() const
   {
      this->mpImpl->applyFft();
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

   void ChebyshevIntegrator::output(MatrixZ& rOut, const bool useReal) const
   {
      this->mpImpl->output(rOut, useReal);
   }

   void ChebyshevIntegrator::outputSpectral(MatrixZ& rOut, const bool useReal) const
   {
      this->mpImpl->outputSpectral(rOut, useReal);
   }

}
}
}
}
