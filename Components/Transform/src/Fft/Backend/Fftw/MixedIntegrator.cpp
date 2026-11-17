/**
 * @file MixedIntegrator.cpp
 * @brief Source of the interface for a generic FFTW based mixed integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Transform/Fft/Backend/Fftw/MixedIntegrator.hpp"
#include "Types/Math.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   void MixedIntegrator::init(const SetupType& setup) const
   {
      //Initialize parent
      IMixedBackend::init(setup);

      int fwdSize = setup.fwdSize();
      int bwdSize = setup.bwdSize();
      int blockSize = setup.blockSize();

      // Set transform scaling
      this->mFftScaling = 1.0/static_cast<MHDFloat>(fwdSize);

      // Get sizes
      this->mBwdSize = bwdSize;
      this->mBlockSize = blockSize;

      // Create the two plans
      const int  *fftSize = &fwdSize;

      // create temporary storage for plan computation
      Matrix    tmpReal(fwdSize, blockSize);
      MatrixZ   tmpCplx(bwdSize, blockSize);

      // Create the real to complex plan
      this->mPlan = fftw_plan_many_dft_r2c(1, fftSize, blockSize, tmpReal.data(), NULL, 1, fwdSize, reinterpret_cast<fftw_complex* >(tmpCplx.data()), NULL, 1, bwdSize, QuICC::Fft::Fftw::Library::planFlag());
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }
   }

   void MixedIntegrator::output(Eigen::Ref<MatrixZ> rOut) const
   {
      rOut *= this->mFftScaling;
   }

   void MixedIntegrator::outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*this->mFftScaling*Math::cI;
         rOut.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*rOut.topRows(this->mSpecSize);
      } else
      {
         MHDFloat sgn = std::pow(-1.0,(order/2)%2)*this->mFftScaling;
         rOut.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*rOut.topRows(this->mSpecSize);
      }
   }

   void MixedIntegrator::outputDiff(Eigen::Ref<MatrixZ> rOut, const int order, const MHDFloat scale, const std::map<int,MHDComplex>& mod) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*this->mFftScaling*Math::cI;
         ArrayZ factor = (sgn*(scale*this->positiveK()).array().pow(order).matrix());
         for(auto m: mod)
         {
            factor(m.first) = m.second*this->mFftScaling;
         }
         rOut.topRows(this->mSpecSize) = factor.asDiagonal()*rOut.topRows(this->mSpecSize);
      } else
      {
         MHDFloat sgn = std::pow(-1.0,(order/2)%2)*this->mFftScaling;
         Array factor = (sgn*(scale*this->positiveK()).array().pow(order).matrix());
         for(auto m: mod)
         {
            assert(m.second.imag() == 0);
            factor(m.first) = m.second.real()*this->mFftScaling;
         }
         rOut.topRows(this->mSpecSize) = factor.asDiagonal()*rOut.topRows(this->mSpecSize);
      }
   }

} // namespace Fftw
} // namespace Backend
} // namespace Fft
} // namespace Transform
} // namespace QuICC
