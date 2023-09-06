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

   void ChebyshevProjector::input(Matrix& tmp, const Matrix& in) const
   {
      this->mpImpl->input(tmp, in);
   }

   void ChebyshevProjector::input(Matrix& tmp, const Matrix& in,
      const int shift) const
   {
      this->mpImpl->input(tmp, in, shift);
   }

   void ChebyshevProjector::input(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const
   {
      this->mpImpl->input(tmp, in, useReal);
   }

   void ChebyshevProjector::input(Matrix& tmp, const MatrixZ& in,
      const int shift, const bool useReal) const
   {
      this->mpImpl->input(tmp, in, shift, useReal);
   }

   void ChebyshevProjector::output(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mpImpl->output(rOut, tmp, useReal);
   }

   void ChebyshevProjector::outputScale(Matrix& rOut) const
   {
      this->mpImpl->outputScale(rOut);
   }

   void ChebyshevProjector::outputScale(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      this->mpImpl->outputScale(rOut, tmp, useReal);
   }

   void ChebyshevProjector::applyFft(Matrix& phys, const Matrix& mods) const
   {
      this->mpImpl->applyFft(phys, mods);
   }

   void ChebyshevProjector::addSolver(const int extraRows) const
   {
      this->mpImpl->addSolver(extraRows);
   }

   void ChebyshevProjector::getSolution(Matrix& tmp, const int zeroRows, const int extraRows, const bool updateSolver) const
   {
      this->mpImpl->getSolution(tmp, zeroRows, extraRows, updateSolver);
   }

   Fftw::DifferentialSolver& ChebyshevProjector::solver() const
   {
      return this->mpImpl->solver();
   }

   Matrix& ChebyshevProjector::getStorage(const StorageKind kind) const
   {
      return this->mpImpl->getStorage(kind);
   }

}
}
}
}
