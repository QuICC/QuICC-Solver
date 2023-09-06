/**
 * @file ChebyshevEnergy.cpp
 * @brief Source of the interface for a generic FFTW based Chebyshev energy reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/ChebyshevEnergy.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   ChebyshevEnergy::ChebyshevEnergy()
   {
   }

   ChebyshevEnergy::~ChebyshevEnergy()
   {
   }

   void ChebyshevEnergy::init(const SetupType& setup) const
   {
      //Initialize parent
      IChebyshevBackend::init(setup);

      int fwdSize = 2*setup.fwdSize();
      int bwdSize = 2*setup.bwdSize();
      int blockSize = setup.blockSize();

      // Set transform scaling
      this->mFftScaling = 1.0/static_cast<MHDFloat>(2*fwdSize);

      // Initialize storage
      this->mTmp.setZero(bwdSize, blockSize);
      this->mTmpComp.setZero(bwdSize, blockSize);
      this->mTmpMid.setZero(fwdSize, blockSize);

      // Set sizes
      this->mFwdSize = setup.fwdSize();
      this->mPadSize = setup.padSize();

      // Compute energy weights
      this->computeEWeights(bwdSize, setup.lower(), setup.upper());

      // Create the two plans
      const int  *fftSize = &fwdSize;

      // Initialise temporary storage
      Matrix tmpF = Matrix::Zero(fwdSize, blockSize);
      Matrix tmpB = Matrix::Zero(bwdSize, blockSize);

      // Create the spectral to physical plan
      const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
      this->mPlan = fftw_plan_many_r2r(1, fftSize, blockSize, tmpB.data(), NULL, 1, bwdSize, tmpF.data(), NULL, 1, fwdSize, bwdKind, QuICC::Fft::Fftw::Library::planFlag());
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }

      // Create the physical to spectral plan
      const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
      this->mFwdPlan = fftw_plan_many_r2r(1, fftSize, blockSize, tmpF.data(), NULL, 1, fwdSize, tmpB.data(), NULL, 1, bwdSize, fwdKind, QuICC::Fft::Fftw::Library::planFlag());
      if(this->mFwdPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }
   }

   void ChebyshevEnergy::computeEWeights(const int size, const MHDFloat lower, const MHDFloat upper) const
   {
      // Initialize energy weights
      this->mEWeights = Array::Zero(size);
      MHDFloat a = (upper - lower)/2.0;
      for(int i = 0; i < size/2; i++)
      {
         MHDFloat n = 2.0*i;
         this->mEWeights(2*i) = 2.0*a*(2.0/(1.0 - n*n));
      }
      this->mEWeights(0) *= 0.5;
   }

   void ChebyshevEnergy::applyPadding(Matrix& rData, const int extraRows) const
   {
      // Set the padded values to zero
      rData.bottomRows(this->mFwdSize+this->mPadSize-extraRows).setZero();
   }

   void ChebyshevEnergy::applyFwdFft(Matrix& mods, const Matrix& phys) const
   {
      fftw_execute_r2r(this->mFwdPlan, const_cast<MHDFloat *>(phys.data()), mods.data());
   }

   void ChebyshevEnergy::setSpectralOperator(const SparseMatrix& mat) const
   {
      this->mSpecOp = mat;
   }

   void ChebyshevEnergy::outputSpectral(Matrix& rOut, const Matrix& tmp) const
   {
      assert(this->mSpecOp.cols() <= tmp.rows());
      assert(this->mSpecOp.rows() <= this->mEWeights.rows());
      rOut.transpose() = this->mFftScaling*this->mEWeights.topRows(this->mSpecOp.rows()).transpose()*this->mSpecOp*tmp.topRows(2*this->mSpecSize);
   }

   void ChebyshevEnergy::output(Matrix& rOut, const Matrix& tmp) const
   {
      int rows = 2*this->mSpecSize;
      rOut.transpose() = this->mFftScaling*this->mEWeights.topRows(rows).transpose()*tmp.topRows(rows);
   }

   void ChebyshevEnergy::square(Matrix& tmp, const Matrix& in,const bool isFirst) const
   {
      if(isFirst)
      {
         tmp = in.array().pow(2);
      } else
      {
         tmp.array() += in.array().pow(2);
      }
   }

   void ChebyshevEnergy::addSolver(const int extraRows) const
   {
      this->mspSolver = std::make_shared<DifferentialSolver>(this->mSpecSize, this->mBlockSize, extraRows);
   }

   void ChebyshevEnergy::getSolution(Matrix& tmp, const int zeroRows, const int extraRows) const
   {
      this->solver().solve(tmp, zeroRows);
      this->applyPadding(tmp, extraRows);
   }

   Matrix& ChebyshevEnergy::getStorage(const StorageKind kind) const
   {
      switch(kind)
      {
         case StorageKind::in:
            return this->mTmp;

         case StorageKind::out:
            return this->mTmpComp;

         default:
            return this->mTmpMid;
      }
   }

   DifferentialSolver& ChebyshevEnergy::solver() const
   {
      return *this->mspSolver;
   }

}
}
}
}
}
