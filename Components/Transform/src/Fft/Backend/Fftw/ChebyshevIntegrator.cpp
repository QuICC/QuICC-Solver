/**
 * @file ChebyshevIntegrator.cpp
 * @brief Source of the interface for a generic FFTW based Chebyshev integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/ChebyshevIntegrator.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   ChebyshevIntegrator::ChebyshevIntegrator()
   {
   }

   ChebyshevIntegrator::~ChebyshevIntegrator()
   {
   }

   void ChebyshevIntegrator::init(const SetupType& setup) const
   {
      //Initialize parent
      IChebyshevBackend::init(setup);

      int fwdSize = setup.fwdSize();
      int bwdSize = setup.bwdSize();
      int blockSize = setup.blockSize();

      // Set transform scaling
      this->mFftScaling = 1.0/static_cast<MHDFloat>(2*fwdSize);

      // Initialize temporary storage
      this->mTmp.setZero(fwdSize, blockSize);
      this->mTmpComp.setZero(fwdSize, blockSize);

      // Set sizes
      this->mBwdSize = setup.bwdSize();

      // Create the plan
      const int  *fftSize = &fwdSize;

      // Initialise temporary storage
      Matrix tmpF = Matrix::Zero(fwdSize, blockSize);
      Matrix tmpB = Matrix::Zero(bwdSize, blockSize);

      // Create the physical to spectral plan
      const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
      this->mPlan = fftw_plan_many_r2r(1, fftSize, blockSize, tmpF.data(), NULL, 1, fwdSize, tmpB.data(), NULL, 1, bwdSize, fwdKind, QuICC::Fft::Fftw::Library::planFlag());
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }
   }

   void ChebyshevIntegrator::setSpectralOperator(const SparseMatrix& mat) const
   {
      this->mSpecOp = mat;
   }

   void ChebyshevIntegrator::setMeanOperator(const SparseMatrix& mat) const
   {
      this->mMeanOp = mat;
   }

   void ChebyshevIntegrator::output(Matrix& rOut) const
   {
      rOut.topRows(this->mSpecSize) = this->mFftScaling*rOut.topRows(this->mSpecSize);
   }

   void ChebyshevIntegrator::outputSpectral(Matrix& rOut) const
   {
      if(this->mMeanOp.size() > 0)
      {
         rOut.block(0, 0, this->mSpecSize, 1) = this->mFftScaling*this->mMeanOp*rOut.block(0,0,this->mMeanOp.cols(),1);
         rOut.block(0, 1, this->mSpecSize, rOut.cols()-1) = this->mFftScaling*this->mSpecOp*rOut.block(0,1,this->mSpecOp.cols(), rOut.cols()-1);
      } else
      {
         rOut.topRows(this->mSpecSize) = this->mFftScaling*this->mSpecOp*rOut.topRows(this->mSpecOp.cols());
      }
   }

   void ChebyshevIntegrator::output(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      if(this->mMeanOp.size() > 0)
      {
         if(useReal)
         {
            rOut.block(0, 0, this->mSpecSize, 1).real() = this->mFftScaling*this->mMeanOp*tmp.block(0,0,this->mMeanOp.cols(),1);
            rOut.block(0, 1, this->mSpecSize, rOut.cols()-1).real() = this->mFftScaling*this->mSpecOp*tmp.block(0,1,this->mSpecOp.cols(), rOut.cols()-1);
         } else
         {
            rOut.block(0, 0, this->mSpecSize, 1).imag() = this->mFftScaling*this->mMeanOp*tmp.block(0,0,this->mMeanOp.cols(),1);
            rOut.block(0, 1, this->mSpecSize, rOut.cols()-1).imag() = this->mFftScaling*this->mSpecOp*tmp.block(0,1,this->mSpecOp.cols(), rOut.cols()-1);
         }
      } else
      {
         if(useReal)
         {
            rOut.topRows(this->mSpecSize).real() = this->mFftScaling*tmp.topRows(this->mSpecSize);
         } else
         {
            rOut.topRows(this->mSpecSize).imag() = this->mFftScaling*tmp.topRows(this->mSpecSize);
         }
      }
   }

   void ChebyshevIntegrator::outputSpectral(MatrixZ& rOut, const Matrix& tmp, const bool useReal) const
   {
      if(useReal)
      {
         rOut.topRows(this->mSpecSize).real() = this->mFftScaling*this->mSpecOp*tmp.topRows(this->mSpecOp.cols());
      } else
      {
         rOut.topRows(this->mSpecSize).imag() = this->mFftScaling*this->mSpecOp*tmp.topRows(this->mSpecOp.cols());
      }
   }

}
}
}
}
}
