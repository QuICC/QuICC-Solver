/**
 * @file ChebyshevProjector.cpp
 * @brief Source of the interface for a generic API for Chebyshev projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/ChebyshevProjector.hpp"

// Project includes
//
#if defined QUICC_FFT_CHEBYSHEV_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/ChebyshevProjector.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_CHEBYSHEV_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/ChebyshevProjector.hpp"
   #define BACKENDIMPL CuFft
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct ChebyshevProjector::BackendImpl: public BACKENDIMPL::ChebyshevProjector
   {
      BackendImpl() = default;
   };

   ChebyshevProjector::ChebyshevProjector()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   ChebyshevProjector::~ChebyshevProjector()
   {
   }

   void ChebyshevProjector::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void ChebyshevProjector::setScaler(const Array& scaler) const
   {
      this->mpImpl->setScaler(scaler);
   }

   void ChebyshevProjector::input(const Matrix& in, const bool needPadding) const
   {
      this->mpImpl->input(in, needPadding);
   }

   void ChebyshevProjector::input(const MatrixZ& in, const bool useReal, const bool needPadding) const
   {
      this->mpImpl->input(in, useReal, needPadding);
   }

   void ChebyshevProjector::io() const
   {
      this->mpImpl->io();
   }

   void ChebyshevProjector::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpImpl->io(out, in);
   }

   void ChebyshevProjector::output(Matrix& rOut) const
   {
      this->mpImpl->output(rOut);
   }

   void ChebyshevProjector::output(MatrixZ& rOut, const bool useReal) const
   {
      this->mpImpl->output(rOut, useReal);
   }

   void ChebyshevProjector::outputScale(Matrix& rOut) const
   {
      this->mpImpl->outputScale(rOut);
   }

   void ChebyshevProjector::outputScale(MatrixZ& rOut, const bool useReal) const
   {
      this->mpImpl->outputScale(rOut, useReal);
   }

   void ChebyshevProjector::applyFft() const
   {
      this->mpImpl->applyFft();
   }

   void ChebyshevProjector::addSolver(const int extraRows) const
   {
      this->mpImpl->addSolver(extraRows);
   }

   void ChebyshevProjector::getSolution(const int zeroRows, const int extraRows, const bool updateSolver) const
   {
      this->mpImpl->getSolution(zeroRows, extraRows, updateSolver);
   }

   Fftw::DifferentialSolver& ChebyshevProjector::solver() const
   {
      return this->mpImpl->solver();
   }

}
}
}
}
