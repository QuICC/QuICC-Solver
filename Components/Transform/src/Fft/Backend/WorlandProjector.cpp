/**
 * @file WorlandProjector.cpp
 * @brief Source of the interface for a generic API for Worland projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/WorlandProjector.hpp"

// Project includes
//
#if defined QUICC_FFT_WORLAND_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/WorlandProjector.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_WORLAND_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/WorlandProjector.hpp"
   #define BACKENDIMPL CuFft
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct WorlandProjector::BackendImpl: public BACKENDIMPL::WorlandProjector
   {
      BackendImpl() = default;
   };

   WorlandProjector::WorlandProjector()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   WorlandProjector::~WorlandProjector()
   {
   }

   void WorlandProjector::init(const SetupType& setup, const int lshift, const int extraN, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      this->mpImpl->init(setup, lshift, extraN, lshiftOnlyParity, alwaysZeroNegative);
   }

   void WorlandProjector::setZFilter(const std::set<int>& filter) const
   {
      this->mpImpl->setZFilter(filter);
   }

   void WorlandProjector::addStorage(const int inExtras, const int outExtras) const
   {
      this->mpImpl->addStorage(inExtras, outExtras);
   }

   void WorlandProjector::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpImpl->io(out, in);
   }

   void WorlandProjector::io(const bool isEven) const
   {
      this->mpImpl->io(isEven);
   }

   void WorlandProjector::input(const Matrix& in, const bool isEven, const bool needPadding) const
   {
      this->mpImpl->input(in, isEven, needPadding);
   }

   void WorlandProjector::input(const MatrixZ& in, const bool isEven, const bool useReal, const bool needPadding) const
   {
      this->mpImpl->input(in, isEven, useReal, needPadding);
   }

   void WorlandProjector::output(Matrix& rOut, const bool isEven) const
   {
      this->mpImpl->output(rOut, isEven);
   }

   void WorlandProjector::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mpImpl->output(rOut, isEven, useReal);
   }

   void WorlandProjector::applyFft() const
   {
      this->mpImpl->applyFft();
   }

   void WorlandProjector::backwardWorland(const bool isEven, const unsigned int id) const
   {
      this->mpImpl->backwardWorland(isEven, id);
   }

   void WorlandProjector::lowerAlpha(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerAlpha(alpha, isEven, id, norm);
   }

   void WorlandProjector::lowerBeta(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerBeta(alpha, isEven, id, norm);
   }

   void WorlandProjector::lowerR2Beta(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerR2Beta(alpha, isEven, id, norm);
   }

   void WorlandProjector::scaleC(const MHDFloat c, const bool isEven, const unsigned int id) const
   {
      this->mpImpl->scaleC(c, isEven, id);
   }

   void WorlandProjector::scaleALPY(const MHDFloat a, const MHDFloat y, const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleALPY(a, y, isEven, lshift, id);
   }

   void WorlandProjector::scaleD(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleD(isEven, lshift, id);
   }

   void WorlandProjector::scaleSphLaplA(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleSphLaplA(isEven, lshift, id);
   }

   void WorlandProjector::scaleSphLaplB(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleSphLaplB(isEven, lshift, id);
   }

   void WorlandProjector::lshift(const unsigned int id, const int lshift, const bool isEven) const
   {
      this->mpImpl->lshift(id, lshift, isEven);
   }

   void WorlandProjector::nshift(const unsigned int id, const int nshift, const bool isEven) const
   {
      this->mpImpl->nshift(id, nshift, isEven);
   }

   void WorlandProjector::copy(const int to, const int from, const int nshift, const bool isEven) const
   {
      this->mpImpl->copy(to, from, nshift, isEven);
   }

   void WorlandProjector::add(const int to, const int from, const int nshift, const bool isEven) const
   {
      this->mpImpl->add(to, from, nshift, isEven);
   }

}
}
}
}
