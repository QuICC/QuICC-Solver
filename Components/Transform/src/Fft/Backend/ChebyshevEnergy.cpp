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

   void ChebyshevEnergy::input(const Matrix& in, const bool needPadding) const
   {
      this->mpImpl->input(in, needPadding);
   }

   void ChebyshevEnergy::input(const MatrixZ& in, const bool useReal, const bool needPadding) const
   {
      this->mpImpl->input(in, useReal, needPadding);
   }

   void ChebyshevEnergy::square(const bool isFirst) const
   {
      this->mpImpl->square(isFirst);
   }

   void ChebyshevEnergy::output(Matrix& rOut) const
   {
      this->mpImpl->output(rOut);
   }

   void ChebyshevEnergy::outputSpectral(Matrix& rOut) const
   {
      this->mpImpl->outputSpectral(rOut);
   }

   void ChebyshevEnergy::io() const
   {
      this->mpImpl->io();
   }

   void ChebyshevEnergy::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpImpl->io(out, in);
   }

   void ChebyshevEnergy::applyFft() const
   {
      this->mpImpl->applyFft();
   }

   void ChebyshevEnergy::applyFwdFft() const
   {
      this->mpImpl->applyFwdFft();
   }

   void ChebyshevEnergy::addSolver(const int extraRows) const
   {
      this->mpImpl->addSolver(extraRows);
   }

   void ChebyshevEnergy::getSolution(const int zeroRows, const int extraRows) const
   {
      this->mpImpl->getSolution(zeroRows, extraRows);
   }

   Fftw::DifferentialSolver& ChebyshevEnergy::solver() const
   {
      return this->mpImpl->solver();
   }

}
}
}
}
