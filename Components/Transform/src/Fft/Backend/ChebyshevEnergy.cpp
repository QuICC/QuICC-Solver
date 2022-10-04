/**
 * @file ChebyshevEnergy.cpp
 * @brief Source of the interface for a generic API for Chebyshev energy reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/ChebyshevEnergy.hpp"

// Project includes
//
#if defined QUICC_FFT_CHEBYSHEV_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/ChebyshevEnergy.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_CHEBYSHEV_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/ChebyshevEnergy.hpp"
   #define BACKENDIMPL CuFft
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct ChebyshevEnergy::BackendImpl: public BACKENDIMPL::ChebyshevEnergy
   {
      BackendImpl() = default;
   };

   ChebyshevEnergy::ChebyshevEnergy()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   ChebyshevEnergy::~ChebyshevEnergy()
   {
   }

   void ChebyshevEnergy::init(const SetupType& setup) const
   {
      this->mpImpl->init(setup);
   }

   void ChebyshevEnergy::setSpectralOperator(const SparseMatrix& mat) const
   {
      this->mpImpl->setSpectralOperator(mat);
   }

   void ChebyshevEnergy::input(Matrix& tmp, const Matrix& in) const
   {
      this->mpImpl->input(tmp, in);
   }

   void ChebyshevEnergy::input(Matrix& tmp, const Matrix& in,
      const int shift) const
   {
      this->mpImpl->input(tmp, in, shift);
   }

   void ChebyshevEnergy::input(Matrix& tmp, const MatrixZ& in,
      const bool useReal) const
   {
      this->mpImpl->input(tmp, in, useReal);
   }

   void ChebyshevEnergy::input(Matrix& tmp, const MatrixZ& in,
      const int shift, const bool useReal) const
   {
      this->mpImpl->input(tmp, in, shift, useReal);
   }

   void ChebyshevEnergy::square(Matrix& tmp, const Matrix& in,const bool isFirst) const
   {
      this->mpImpl->square(tmp, in, isFirst);
   }

   void ChebyshevEnergy::output(Matrix& rOut, const Matrix& tmp) const
   {
      this->mpImpl->output(rOut, tmp);
   }

   void ChebyshevEnergy::outputSpectral(Matrix& rOut, const Matrix& tmp) const
   {
      this->mpImpl->outputSpectral(rOut, tmp);
   }

   void ChebyshevEnergy::applyFft(Matrix& phys, const Matrix& mods) const
   {
      this->mpImpl->applyFft(phys, mods);
   }

   void ChebyshevEnergy::applyFwdFft(Matrix& mods, const Matrix& phys) const
   {
      this->mpImpl->applyFwdFft(mods, phys);
   }

   void ChebyshevEnergy::addSolver(const int extraRows) const
   {
      this->mpImpl->addSolver(extraRows);
   }

   void ChebyshevEnergy::getSolution(Matrix& tmp, const int zeroRows, const int extraRows) const
   {
      this->mpImpl->getSolution(tmp, zeroRows, extraRows);
   }

   Matrix& ChebyshevEnergy::getStorage(const StorageKind kind) const
   {
      return this->mpImpl->getStorage(kind);
   }

   Fftw::DifferentialSolver& ChebyshevEnergy::solver() const
   {
      return this->mpImpl->solver();
   }

}
}
}
}
