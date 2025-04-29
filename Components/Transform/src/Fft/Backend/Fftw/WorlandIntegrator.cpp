/**
 * @file WorlandIntegrator.cpp
 * @brief Source of the interface for a generic FFTW based Worland integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/WorlandIntegrator.hpp"

// Project includes
//
#include "Types/Math.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   WorlandIntegrator::WorlandIntegrator()
   {
   }

   WorlandIntegrator::~WorlandIntegrator()
   {
   }

   void WorlandIntegrator::init(const SetupType& setup, const int lshift, const int extraN, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      //Initialize parent
      IWorlandBackend::init(setup, lshift, extraN, lshiftOnlyParity, alwaysZeroNegative);
      this->mpWorkTmp = &this->mOutTmp;

      int fwdSize = setup.fwdSize();
      int bwdSize = setup.bwdSize();

      // Set internal spectral resolution
      this->setWSize();
      assert(this->mWSize <= bwdSize);

      // Set transform scaling
      this->mFftScaling = std::sqrt(Math::PI)/static_cast<MHDFloat>(2*fwdSize);

      // Initialize main temporary storage
      int blockSize = std::max(this->mEBlockSize,this->mOBlockSize);
      // Input temporary storage
      this->initStorage(fwdSize, blockSize, 1, this->mInTmp);
      // Output temporary storage
      this->initStorage(bwdSize, blockSize, 1, this->mOutTmp);

      // Set sizes
      this->mBwdSize = setup.bwdSize();

      // Create the plan size
      const int  *fftSize = &fwdSize;

      // Initialise temporary even storage
      blockSize = this->mEBlockSize;
      Matrix tmpF = Matrix::Zero(fwdSize, blockSize);
      Matrix tmpB = Matrix::Zero(bwdSize, blockSize);

      // Create the even physical to spectral plan
      const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
      this->mPlan = fftw_plan_many_r2r(1, fftSize, blockSize, tmpF.data(), NULL, 1, fwdSize, tmpB.data(), NULL, 1, bwdSize, fwdKind, QuICC::Fft::Fftw::Library::planFlag());
      if(this->mPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }

      // Initialise temporary odd storage
      blockSize = this->mOBlockSize;
      tmpF = Matrix::Zero(fwdSize, blockSize);
      tmpB = Matrix::Zero(bwdSize, blockSize);

      // Create the odd physical to spectral plan
      const fftw_r2r_kind ofwdKind[] = {FFTW_REDFT11};
      this->mOddPlan = fftw_plan_many_r2r(1, fftSize, blockSize, tmpF.data(), NULL, 1, fwdSize, tmpB.data(), NULL, 1, bwdSize, ofwdKind, QuICC::Fft::Fftw::Library::planFlag());
      if(this->mOddPlan == NULL)
      {
         throw  std::logic_error("FFTW plan failed!");
      }

      // Initialize Jacobi shift matrices
      this->initJ();

      // Initialize Banded matrices
      this->initBanded(false);
      this->initBanded(true);
   }

   int WorlandIntegrator::lSize(const int l) const
   {
      return this->mWSize - l/2;
   }

   int WorlandIntegrator::lSize(const int i, const int l) const
   {
      assert(l>=0);
      return this->mWSize - i;
   }

   void WorlandIntegrator::io(const bool isEven) const
   {
      this->setPlan(isEven);

      this->io(this->mOutTmp.at(0).data(), this->mInTmp.at(0).data());
   }

   void WorlandIntegrator::input(const Matrix& in, const bool isEven) const
   {
      Matrix& inTmp = this->mInTmp.at(0);
      int start = 0;
      for(auto loc: *this->pLoc(isEven))
      {
         inTmp.block(0, start, inTmp.rows(), std::get<2>(loc)) = in.block(0, std::get<1>(loc), in.rows(), std::get<2>(loc));
         start += std::get<2>(loc);
      }
   }

   void WorlandIntegrator::input(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      Matrix& inTmp = this->mInTmp.at(0);
      int start = 0;

      if(useReal)
      {
         for(auto loc: *this->pLoc(isEven))
         {
            inTmp.block(0, start, inTmp.rows(), std::get<2>(loc)) = in.block(0, std::get<1>(loc), in.rows(), std::get<2>(loc)).real();
            start += std::get<2>(loc);
         }
      } else
      {
         for(auto loc: *this->pLoc(isEven))
         {
            inTmp.block(0, start, inTmp.rows(), std::get<2>(loc)) = in.block(0, std::get<1>(loc), in.rows(), std::get<2>(loc)).imag();
            start += std::get<2>(loc);
         }
      }
   }

   void WorlandIntegrator::applyFft() const
   {
      fftw_execute_r2r(*this->mpPlan, const_cast<MHDFloat *>(this->mpIn), this->mpOut);
   }

   void WorlandIntegrator::partialForwardWorland(const int l, const int i0, const int start, const int cols, const std::vector<Matrix>& banded, Matrix& out) const
   {
      for(int i = i0; i < l/2; ++i)
      {
         this->applyPair(out, start, cols, this->lSize(i,l), banded.at(i));
      }
   }

   void WorlandIntegrator::forwardWorland(const bool isEven, const int id) const
   {
      // Reset current l
      this->resetLocations(isEven, id);

      Matrix& outTmp = this->mOutTmp.at(id);
      const std::vector<Matrix>& banded = this->banded(isEven);

      int start = 0;
      int i0 = 0;
      int cols = outTmp.cols();
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         cols = outTmp.cols() - start;
         if(l < 2)
         {
            outTmp.block(0, start, this->lSize(l), std::get<2>(loc)).array() *= this->mFftScaling;
            if(l == 0)
            {
               outTmp.block(0, start, 1, std::get<2>(loc)).array() *= 1.0/std::sqrt(2.0);
            }
         } else
         {
            if(i0 == 0)
            {
               outTmp.block(0, start, this->lSize(0,l), cols) *= this->mFftScaling;
            }
            this->partialForwardWorland(l, i0, start, cols, banded, outTmp);
            i0 = l/2;
         }
         std::get<3>(loc) = this->lSize(l);
         start += std::get<2>(loc);
      }
   }

   void WorlandIntegrator::lowerBeta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      Matrix& outTmp = this->mOutTmp.at(id);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         Matrix V;
         this->buildShiftV(V, this->lSize(l), alpha, l-0.5, norm);
         this->applyTriSolve(outTmp, start, cols, V);
         if(l == 1)
         {
            outTmp.block(0, start, 1, std::get<2>(loc)).array() *= 1.0/std::sqrt(2.0);
         }
         std::get<4>(loc)--;
         start += cols;
      }
   }

   void WorlandIntegrator::raiseBeta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      Matrix& outTmp = this->mOutTmp.at(id);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l < 0)
         {
            outTmp.block(0, start, outTmp.rows(), cols).setZero();
         } else
         {
            Matrix V;
            this->buildShiftV(V, this->lSize(l), alpha, l+0.5, norm);
            if(l==0)
            {
               V(1,0) *= std::sqrt(2.0);
            }
            this->applyTriProduct(outTmp, start, cols, V);
            std::get<4>(loc)++;
         }
         std::get<3>(loc)--;
         start += cols;
      }
   }

   void WorlandIntegrator::lowerR2Beta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      Matrix& outTmp = this->mOutTmp.at(id);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l < 0)
         {
            outTmp.block(0, start, outTmp.rows(), cols).setZero();
         } else
         {
            Matrix M;
            this->buildShiftM(M, this->lSize(l), alpha, l-0.5, norm);
            this->applyTriProduct(outTmp, start, cols, M);
            std::get<4>(loc)--;
         }
         start += cols;
      }
   }

   void WorlandIntegrator::raiseR2Beta(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm, const bool scaleL0) const
   {
      Matrix& outTmp = this->mOutTmp.at(id);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         Matrix M;
         this->buildShiftM(M, this->lSize(l), alpha, l+0.5, norm);
         this->applyTriSolve(outTmp, start, cols, M);
         if(scaleL0 && l == 0)
         {
            outTmp.block(0, start, 1, std::get<2>(loc)).array() *= std::sqrt(2.0);
         }
         std::get<4>(loc)++;
         start += cols;
      }
   }

   void WorlandIntegrator::lowerAlpha(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      Matrix& outTmp = this->mOutTmp.at(id);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         Matrix U;
         this->buildShiftU(U, this->lSize(l), alpha, l-0.5, norm);
         this->applyTriSolve(outTmp, start, cols, U);
         if(l == 1)
         {
            outTmp.block(0, start, 1, std::get<2>(loc)).array() *= 1.0/std::sqrt(2.0);
         }
         start += cols;
      }
   }

   void WorlandIntegrator::raiseAlpha(const MHDFloat alpha, const bool isEven, const int id, const MHDFloat norm) const
   {
      throw std::logic_error("Raise alpha operator has not been tested!");

      Matrix& outTmp = this->mOutTmp.at(id);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l < 0)
         {
            outTmp.block(0, start, outTmp.rows(), cols).setZero();
         } else
         {
            Matrix U;
            this->buildShiftU(U, this->lSize(l), alpha, l-0.5, norm);
            this->applyTriProduct(outTmp, start, cols, U);
         }
         std::get<3>(loc)--;
         start += cols;
      }
   }

   void WorlandIntegrator::applyI2(const bool isEven, const int id) const
   {
      Matrix& outTmp = this->mOutTmp.at(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         int r = this->lSize(l);
         Internal::MHDFloat a = -MHD_MP(0.5);
         Internal::MHDFloat b = -MHD_MP(0.5);
         ::QuICC::SparseSM::Worland::I2 spasm(r, r, a, b, l);

         Matrix bd;
         spasm.buildOp(bd);
         this->applyBandProduct(outTmp, start, cols, bd);

         start += cols;
      }
   }

   void WorlandIntegrator::applyI4(const bool isEven, const int id) const
   {
      Matrix& outTmp = this->mOutTmp.at(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         int r = this->lSize(l);
         Internal::MHDFloat a = -MHD_MP(0.5);
         Internal::MHDFloat b = -MHD_MP(0.5);
         ::QuICC::SparseSM::Worland::I4 spasm(r, r, a, b, l);

         Matrix bd;
         spasm.buildOp(bd);
         this->applyBandProduct(outTmp, start, cols, bd);

         start += cols;
      }
   }

   void WorlandIntegrator::output(Matrix& rOut, const bool isEven) const
   {
      Matrix& outTmp = this->mOutTmp.at(0);
      int start = 0;
      for(auto loc: *this->pLoc(isEven))
      {
         int cols = std::get<2>(loc);
         rOut.block(0, std::get<1>(loc), this->mSpecSize, cols) = outTmp.block(0, start, this->mSpecSize, cols);
         start += cols;
      }

      // Zero unused values
      for(auto loc: this->mZLoc)
      {
         int s = std::get<1>(loc);
         int cols = std::get<2>(loc);
         rOut.block(0, s, this->mSpecSize, cols).setZero();
      }
   }

   void WorlandIntegrator::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      Matrix& outTmp = this->mOutTmp.at(0);
      int start = 0;
      if(useReal)
      {
         for(auto loc: *this->pLoc(isEven))
         {
            int cols = std::get<2>(loc);
            rOut.block(0, std::get<1>(loc), this->mSpecSize, cols).real() = outTmp.block(0, start, this->mSpecSize, cols);
            start += cols;
         }

         // Zero unused values
         for(auto loc: this->mZLoc)
         {
            int s = std::get<1>(loc);
            int cols = std::get<2>(loc);
            rOut.block(0, s, this->mSpecSize, cols).real().setZero();
         }
      } else
      {
         for(auto loc: *this->pLoc(isEven))
         {
            int cols = std::get<2>(loc);
            rOut.block(0, std::get<1>(loc), this->mSpecSize, cols).imag() = outTmp.block(0, start, this->mSpecSize, cols);
            start += cols;
         }

         // Zero unused values
         for(auto loc: this->mZLoc)
         {
            int s = std::get<1>(loc);
            int cols = std::get<2>(loc);
            rOut.block(0, s, this->mSpecSize, cols).imag().setZero();
         }
      }

   }

   void WorlandIntegrator::initBanded(const bool isEven) const
   {
      if(this->pLoc(isEven,0)->size() > 0)
      {
         const Matrix& J = this->J(isEven);

         const MHDFloat normV = 1.0;
         const MHDFloat normM = 1.0;
         std::vector<Matrix>* pBanded;
         if(isPhysEven(isEven))
         {
            pBanded = &this->mEBanded;
         } else
         {
            pBanded = &this->mOBanded;
         }
         int l = std::get<4>(this->pLoc(isEven,0)->back());
         if(l >= 2)
         {
            for(int i = 0; i < l/2; ++i)
            {
               Matrix PS;
               this->buildShiftPair(PS, true, i, this->lSize(i,l), J, normV, normM);
               pBanded->push_back(PS);
            }
         }
      }
   }

}
}
}
}
}
