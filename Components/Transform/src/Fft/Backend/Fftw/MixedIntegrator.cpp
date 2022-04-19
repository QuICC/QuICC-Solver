/**
 * @file MixedIntegrator.cpp
 * @brief Source of the interface for a generic FFTW based mixed integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/MixedIntegrator.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   MixedIntegrator::MixedIntegrator()
      : mOutMap(NULL,0,0)
   {
   }

   MixedIntegrator::~MixedIntegrator()
   {
   }

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
      this->mPlan = fftw_plan_many_dft_r2c(1, fftSize, blockSize, tmpReal.data(), NULL, 1, fwdSize, reinterpret_cast<fftw_complex* >(tmpCplx.data()), NULL, 1, bwdSize, Library::planFlag());
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }
   }

   void MixedIntegrator::output(MatrixZ& rOut) const
   {
      rOut.topRows(this->mSpecSize) = this->mFftScaling*this->mOutMap.topRows(this->mSpecSize);
   }

   void MixedIntegrator::outputDiff(MatrixZ& rOut, const int order, const MHDFloat scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*this->mFftScaling*Math::cI;
         rOut.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*this->mOutMap.topRows(this->mSpecSize);
      } else
      {
         MHDFloat sgn = std::pow(-1.0,(order/2)%2)*this->mFftScaling;
         rOut.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*this->mOutMap.topRows(this->mSpecSize);
      }
   }

   void MixedIntegrator::outputDiff(MatrixZ& rOut, const int order, const MHDFloat scale, const std::map<int,MHDComplex>& mod) const
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
         rOut.topRows(this->mSpecSize) = factor.asDiagonal()*this->mOutMap.topRows(this->mSpecSize);
      } else
      {
         MHDFloat sgn = std::pow(-1.0,(order/2)%2)*this->mFftScaling;
         Array factor = (sgn*(scale*this->positiveK()).array().pow(order).matrix());
         for(auto m: mod)
         {
            assert(m.second.imag() == 0);
            factor(m.first) = m.second.real()*this->mFftScaling;
         }
         rOut.topRows(this->mSpecSize) = factor.asDiagonal()*this->mOutMap.topRows(this->mSpecSize);
      }
   }

   void MixedIntegrator::io(MatrixZ& rOut, const Matrix& in) const
   {
      if(rOut.rows() == this->mBwdSize)
      {
         this->io(rOut.data(), in.data());
      } else
      {
         this->mTmp.resize(this->mBwdSize, this->mBlockSize);
         this->io(this->mTmp.data(), in.data());
      }
   }

   void MixedIntegrator::io(MHDComplex* out, const MHDFloat* in) const
   {
      this->mpOut = out;
      this->mpIn = in;

      new (&this->mOutMap) Eigen::Map<MatrixZ>(this->mpOut, this->mBwdSize, this->mBlockSize);
   }

   void MixedIntegrator::applyFft() const
   {
      fftw_execute_dft_r2c(this->mPlan, const_cast<MHDFloat*>(this->mpIn), reinterpret_cast<fftw_complex* >(this->mpOut));
   }

}
}
}
}
}
