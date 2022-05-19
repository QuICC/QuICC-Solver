/**
 * @file WorlandIntegrator.cpp
 * @brief Source of the interface for a generic API for Worland integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/WorlandIntegrator.hpp"

// Project includes
//
#if defined QUICC_FFT_WORLAND_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/WorlandIntegrator.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_WORLAND_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/WorlandIntegrator.hpp"
   #define BACKENDIMPL CuFft
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct WorlandIntegrator::BackendImpl: public BACKENDIMPL::WorlandIntegrator
   {
      BackendImpl() = default;
   };

   WorlandIntegrator::WorlandIntegrator()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   WorlandIntegrator::~WorlandIntegrator()
   {
   }

   void WorlandIntegrator::init(const SetupType& setup, const int lshift, const int extraN, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      this->mpImpl->init(setup, lshift, extraN, lshiftOnlyParity, alwaysZeroNegative);
   }

   void WorlandIntegrator::setZFilter(const std::set<int>& filter) const
   {
      this->mpImpl->setZFilter(filter);
   }

   void WorlandIntegrator::addStorage(const int inExtras, const int outExtras) const
   {
      this->mpImpl->addStorage(inExtras, outExtras);
   }

   void WorlandIntegrator::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpImpl->io(out, in);
   }

   void WorlandIntegrator::io(const bool isEven) const
   {
      this->mpImpl->io(isEven);
   }

   void WorlandIntegrator::input(const Matrix& in, const bool isEven) const
   {
      this->mpImpl->input(in, isEven);
   }

   void WorlandIntegrator::input(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mpImpl->input(in, isEven, useReal);
   }

   void WorlandIntegrator::applyFft() const
   {
      this->mpImpl->applyFft();
   }

   void WorlandIntegrator::forwardWorland(const bool isEven, const int id) const
   {
      this->mpImpl->forwardWorland(isEven, id);
   }

   void WorlandIntegrator::lowerBeta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerBeta(alpha, isEven, id, norm);
   }

   void WorlandIntegrator::raiseBeta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->raiseBeta(alpha, isEven, id, norm);
   }

   void WorlandIntegrator::lowerR2Beta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerR2Beta(alpha, isEven, id, norm);
   }

   void WorlandIntegrator::raiseR2Beta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->raiseR2Beta(alpha, isEven, id, norm);
   }

   void WorlandIntegrator::lowerAlpha(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerAlpha(alpha, isEven, id, norm);
   }

   void WorlandIntegrator::raiseAlpha(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->raiseAlpha(alpha, isEven, id, norm);
   }

   void WorlandIntegrator::applyI2(const bool isEven, const int id) const
   {
      this->mpImpl->applyI2(isEven, id);
   }

   void WorlandIntegrator::applyI4(const bool isEven, const int id) const
   {
      this->mpImpl->applyI4(isEven, id);
   }

   void WorlandIntegrator::output(Matrix& rOut, const bool isEven) const
   {
      this->mpImpl->output(rOut, isEven);
   }

   void WorlandIntegrator::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mpImpl->output(rOut, isEven, useReal);
   }

   void WorlandIntegrator::scaleC(const MHDFloat c, const bool isEven, const unsigned int id) const
   {
      this->mpImpl->scaleC(c, isEven, id);
   }

   void WorlandIntegrator::scaleALPY(const MHDFloat a, const MHDFloat y, const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleALPY(a, y, isEven, lshift, id);
   }

   void WorlandIntegrator::scaleD(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleD(isEven, lshift, id);
   }

   void WorlandIntegrator::scaleSphLaplA(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleSphLaplA(isEven, lshift, id);
   }

   void WorlandIntegrator::scaleSphLaplB(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleSphLaplB(isEven, lshift, id);
   }

   void WorlandIntegrator::lshift(const unsigned int id, const int lshift, const bool isEven) const
   {
      this->mpImpl->lshift(id, lshift, isEven);
   }

   void WorlandIntegrator::nshift(const unsigned int id, const int nshift, const bool isEven) const
   {
      this->mpImpl->nshift(id, nshift, isEven);
   }

   void WorlandIntegrator::copy(const int to, const int from, const int nshift, const bool isEven) const
   {
      this->mpImpl->copy(to, from, nshift, isEven);
   }

   void WorlandIntegrator::add(const int to, const int from, const int nshift, const bool isEven) const
   {
      this->mpImpl->add(to, from, nshift, isEven);
   }

}
}
}
}
