/**
 * @file ALegendreProjector.cpp
 * @brief Source of the interface for a generic API for ALegendre projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/ALegendreProjector.hpp"

// Project includes
//
#if defined QUICC_FFT_ALEGENDRE_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/ALegendreProjector.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_ALEGENDRE_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/ALegendreProjector.hpp"
   #define BACKENDIMPL CuFft
#elif defined QUICC_FFT_ALEGENDRE_PFSOLVE
   #include "QuICC/Transform/Fft/Backend/PfSolve/ALegendreProjector.hpp"
   #define BACKENDIMPL PfSolve_parallALT
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct ALegendreProjector::BackendImpl: public BACKENDIMPL::ALegendreProjector
   {
      BackendImpl() = default;
   };

   ALegendreProjector::ALegendreProjector()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   ALegendreProjector::~ALegendreProjector()
   {
   }

   void ALegendreProjector::init(const SetupType& setup, const int lshift, const int extraN, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      this->mpImpl->init(setup, lshift, extraN, lshiftOnlyParity, alwaysZeroNegative);
   }

   void ALegendreProjector::setZFilter(const std::set<int>& filter) const
   {
      this->mpImpl->setZFilter(filter);
   }

   void ALegendreProjector::addStorage(const int inExtras, const int outExtras) const
   {
      this->mpImpl->addStorage(inExtras, outExtras);
   }

   void ALegendreProjector::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpImpl->io(out, in);
   }

   void ALegendreProjector::io(const bool isEven) const
   {
      this->mpImpl->io(isEven);
   }

   void ALegendreProjector::input(const Matrix& in, const bool isEven, const bool needPadding) const
   {
      this->mpImpl->input(in, isEven, needPadding);
   }

   void ALegendreProjector::input(const MatrixZ& in, const bool isEven, const bool useReal, const bool needPadding) const
   {
      this->mpImpl->input(in, isEven, useReal, needPadding);
   }

   void ALegendreProjector::output(Matrix& rOut, const bool isEven) const
   {
      this->mpImpl->output(rOut, isEven);
   }

   void ALegendreProjector::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mpImpl->output(rOut, isEven, useReal);
   }

   void ALegendreProjector::applyFft() const
   {
      this->mpImpl->applyFft();
   }

   void ALegendreProjector::backwardALegendre(const bool isEven, const unsigned int id) const
   {
      this->mpImpl->backwardALegendre(isEven, id);
   }

   void ALegendreProjector::lowerAlpha(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerAlpha(alpha, isEven, id, norm);
   }

   void ALegendreProjector::lowerBeta(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerBeta(alpha, isEven, id, norm);
   }

   void ALegendreProjector::lowerR2Beta(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerR2Beta(alpha, isEven, id, norm);
   }

   void ALegendreProjector::scaleC(const MHDFloat c, const bool isEven, const unsigned int id) const
   {
      this->mpImpl->scaleC(c, isEven, id);
   }

   void ALegendreProjector::scaleALPY(const MHDFloat a, const MHDFloat y, const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleALPY(a, y, isEven, lshift, id);
   }

   void ALegendreProjector::scaleD(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleD(isEven, lshift, id);
   }

   void ALegendreProjector::scaleSphLaplA(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleSphLaplA(isEven, lshift, id);
   }

   void ALegendreProjector::scaleSphLaplB(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleSphLaplB(isEven, lshift, id);
   }

   void ALegendreProjector::lshift(const unsigned int id, const int lshift, const bool isEven) const
   {
      this->mpImpl->lshift(id, lshift, isEven);
   }

   void ALegendreProjector::nshift(const unsigned int id, const int nshift, const bool isEven) const
   {
      this->mpImpl->nshift(id, nshift, isEven);
   }

   void ALegendreProjector::copy(const int to, const int from, const int nshift, const bool isEven) const
   {
      this->mpImpl->copy(to, from, nshift, isEven);
   }

   void ALegendreProjector::add(const int to, const int from, const int nshift, const bool isEven) const
   {
      this->mpImpl->add(to, from, nshift, isEven);
   }

}
}
}
}
