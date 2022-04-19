/**
 * @file ChebyshevProjector.cpp
 * @brief Source of the interface for a generic FFTW based Chebyshev projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/ChebyshevProjector.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   ChebyshevProjector::ChebyshevProjector()
   {
   }

   ChebyshevProjector::~ChebyshevProjector()
   {
   }

   void ChebyshevProjector::init(const SetupType& setup) const
   {
      //Initialize parent
      IChebyshevBackend::init(setup);

      int fwdSize = setup.fwdSize();
      int bwdSize = setup.bwdSize();
      int blockSize = setup.blockSize();

      // Create the plan
      const int  *fftSize = &fwdSize;

      // Initialize storage
      this->mTmp.setZero(bwdSize, blockSize);
      this->mTmpComp.setZero(bwdSize, blockSize);

      // Set sizes
      this->mPadSize = setup.padSize();

      // Initialise temporary storage
      Matrix tmpF = Matrix::Zero(fwdSize, blockSize);
      Matrix tmpB = Matrix::Zero(bwdSize, blockSize);

      // Create the spectral to physical plan
      const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
      this->mPlan = fftw_plan_many_r2r(1, fftSize, blockSize, tmpB.data(), NULL, 1, bwdSize, tmpF.data(), NULL, 1, fwdSize, bwdKind, Library::planFlag());
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }
   }

   void ChebyshevProjector::input(const Matrix& in, const bool needPadding) const
   {
      this->mTmp.topRows(this->mSpecSize) = in.topRows(this->mSpecSize);

      // Apply padding if required
      if(needPadding)
      {
         this->applyPadding(this->mTmp);
      }
   }

   void ChebyshevProjector::setScaler(const Array& scaler) const
   {
      this->mScaler = scaler;
   }

   void ChebyshevProjector::input(const MatrixZ& in, const bool useReal, const bool needPadding) const
   {
      if(useReal)
      {
         this->mTmp.topRows(this->mSpecSize) = in.real().topRows(this->mSpecSize);
      } else
      {
         this->mTmp.topRows(this->mSpecSize) = in.imag().topRows(this->mSpecSize);
      }

      // Apply padding if required
      if(needPadding)
      {
         this->applyPadding(this->mTmp);
      }
   }

   void ChebyshevProjector::io() const
   {
      this->io(this->mTmpComp.data(), this->mTmp.data());
   }

   void ChebyshevProjector::output(Matrix& rOut) const
   {
      this->io(rOut.data(), this->mTmp.data());
   }

   void ChebyshevProjector::output(MatrixZ& rOut, const bool useReal) const
   {
      if(useReal)
      {
         rOut.real() = this->mTmpComp;
      } else
      {
         rOut.imag() = this->mTmpComp;
      }
   }

   void ChebyshevProjector::outputScale(Matrix& rOut) const
   {
      rOut = this->mScaler.asDiagonal()*rOut;
   }

   void ChebyshevProjector::outputScale(MatrixZ& rOut, const bool useReal) const
   {
      if(useReal)
      {
         rOut.real() = this->mScaler.asDiagonal()*this->mTmpComp;
      } else
      {
         rOut.imag() = this->mScaler.asDiagonal()*this->mTmpComp;
      }
   }

   void ChebyshevProjector::applyPadding(Matrix& rData, const int extraRows) const
   {
      if(extraRows >= 0)
      {
         // Set the padded values to zero
         rData.bottomRows(this->mPadSize-extraRows).setZero();
      } else
      {
         rData.bottomRows(-extraRows).setZero();
      }
   }

   void ChebyshevProjector::applyFft() const
   {
      fftw_execute_r2r(this->mPlan, const_cast<MHDFloat *>(this->mpIn), this->mpOut);
   }

   void ChebyshevProjector::addSolver(const int extraRows) const
   {
      this->mspSolver = std::make_shared<DifferentialSolver>(this->mSpecSize, this->mBlockSize, extraRows);
   }

   void ChebyshevProjector::getSolution(const int zeroRows, const int extraRows, const bool updateSolver) const
   {
      this->solver().solve(zeroRows, this->mTmp);
      this->applyPadding(this->mTmp, extraRows);
      if(updateSolver)
      {
         this->solver().input(this->mTmp, 0);
      }
   }

   DifferentialSolver& ChebyshevProjector::solver() const
   {
      return *this->mspSolver;
   }

}
}
}
}
}
