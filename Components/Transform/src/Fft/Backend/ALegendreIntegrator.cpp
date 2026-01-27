/**
 * @file ALegendreIntegrator.cpp
 * @brief Source of the interface for a generic API for ALegendre integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/ALegendreIntegrator.hpp"

// Project includes
//
#if defined QUICC_FFT_ALEGENDRE_FFTW
   #include "QuICC/Transform/Fft/Backend/Fftw/ALegendreIntegrator.hpp"
   #define BACKENDIMPL Fftw
#elif defined QUICC_FFT_ALEGENDRE_CUFFT
   #include "QuICC/Transform/Fft/Backend/CuFft/ALegendreIntegrator.hpp"
   #define BACKENDIMPL CuFft
#elif defined QUICC_FFT_ALEGENDRE_PFSOLVE
   #include "QuICC/Transform/Fft/Backend/PfSolve/ALegendreIntegrator.hpp"
   #define BACKENDIMPL PfSolve_parallALT
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

   struct ALegendreIntegrator::BackendImpl: public BACKENDIMPL::ALegendreIntegrator
   {
      BackendImpl() = default;
   };

   ALegendreIntegrator::ALegendreIntegrator()
   {
      this->mpImpl = std::make_shared<BackendImpl>();
   }

   ALegendreIntegrator::~ALegendreIntegrator()
   {
   }

   void ALegendreIntegrator::init(const SetupType& setup, const int lshift, const int extraN, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      this->mpImpl->init(setup, lshift, extraN, lshiftOnlyParity, alwaysZeroNegative);
   }

   void ALegendreIntegrator::setZFilter(const std::set<int>& filter) const
   {
      this->mpImpl->setZFilter(filter);
   }

   void ALegendreIntegrator::addStorage(const int inExtras, const int outExtras) const
   {
      this->mpImpl->addStorage(inExtras, outExtras);
   }

   void ALegendreIntegrator::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpImpl->io(out, in);
   }

   void ALegendreIntegrator::io(const bool isEven) const
   {
      this->mpImpl->io(isEven);
   }

   void ALegendreIntegrator::input(const Matrix& in, const bool isEven) const
   {
      this->mpImpl->input(in, isEven);
   }

   void ALegendreIntegrator::input(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      this->mpImpl->input(in, isEven, useReal);
   }

   void ALegendreIntegrator::applyFft() const
   {
      this->mpImpl->applyFft();
   }

   void ALegendreIntegrator::forwardALegendre(const bool isEven, const int id) const
   {
      this->mpImpl->forwardALegendre(isEven, id);
   }

   void ALegendreIntegrator::lowerBeta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerBeta(alpha, isEven, id, norm);
   }

   void ALegendreIntegrator::raiseBeta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->raiseBeta(alpha, isEven, id, norm);
   }

   void ALegendreIntegrator::lowerR2Beta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerR2Beta(alpha, isEven, id, norm);
   }

   void ALegendreIntegrator::raiseR2Beta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm, const bool scaleL0) const
   {
      this->mpImpl->raiseR2Beta(alpha, isEven, id, norm, scaleL0);
   }

   void ALegendreIntegrator::lowerAlpha(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->lowerAlpha(alpha, isEven, id, norm);
   }

   void ALegendreIntegrator::raiseAlpha(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      this->mpImpl->raiseAlpha(alpha, isEven, id, norm);
   }

   void ALegendreIntegrator::applyI2(const bool isEven, const int id) const
   {
      this->mpImpl->applyI2(isEven, id);
   }

   void ALegendreIntegrator::applyI4(const bool isEven, const int id) const
   {
      this->mpImpl->applyI4(isEven, id);
   }

   void ALegendreIntegrator::output(Matrix& rOut, const bool isEven) const
   {
      this->mpImpl->output(rOut, isEven);
   }

   void ALegendreIntegrator::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      this->mpImpl->output(rOut, isEven, useReal);
   }

   void ALegendreIntegrator::scaleC(const MHDFloat c, const bool isEven, const unsigned int id) const
   {
      this->mpImpl->scaleC(c, isEven, id);
   }

   void ALegendreIntegrator::scaleALPY(const MHDFloat a, const MHDFloat y, const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleALPY(a, y, isEven, lshift, id);
   }

   void ALegendreIntegrator::scaleD(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleD(isEven, lshift, id);
   }

   void ALegendreIntegrator::scaleSphLaplA(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleSphLaplA(isEven, lshift, id);
   }

   void ALegendreIntegrator::scaleSphLaplB(const bool isEven, const int lshift, const unsigned int id) const
   {
      this->mpImpl->scaleSphLaplB(isEven, lshift, id);
   }

   void ALegendreIntegrator::lshift(const unsigned int id, const int lshift, const bool isEven) const
   {
      this->mpImpl->lshift(id, lshift, isEven);
   }

   void ALegendreIntegrator::nshift(const unsigned int id, const int nshift, const bool isEven) const
   {
      this->mpImpl->nshift(id, nshift, isEven);
   }

   void ALegendreIntegrator::copy(const int to, const int from, const int nshift, const bool isEven) const
   {
      this->mpImpl->copy(to, from, nshift, isEven);
   }

   void ALegendreIntegrator::add(const int to, const int from, const int nshift, const bool isEven) const
   {
      this->mpImpl->add(to, from, nshift, isEven);
   }

}
}
}
}
