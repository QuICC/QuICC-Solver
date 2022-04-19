/**
 * @file WorlandProjector.cpp
 * @brief Source of the interface for a generic FFTW based Worland projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/WorlandProjector.hpp"

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   WorlandProjector::WorlandProjector()
   {
   }

   WorlandProjector::~WorlandProjector()
   {
   }

   void WorlandProjector::init(const SetupType& setup, const int lshift, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      //Initialize parent
      IWorlandBackend::init(setup, lshift, lshiftOnlyParity, alwaysZeroNegative);
      this->mpWorkTmp = &this->mInTmp;

      int fwdSize = setup.fwdSize();
      int bwdSize = setup.bwdSize();

      // Set internal spectral resolution
      this->setWSize(0);

      // Create the plan
      const int  *fftSize = &fwdSize;

      // Initialize temporary storage
      int blockSize = std::max(this->mEBlockSize,this->mOBlockSize);
      // Input temporary storage
      this->initStorage(bwdSize, blockSize, 1, this->mInTmp);
      // Output temporary storage
      this->initStorage(fwdSize, blockSize, 1, this->mOutTmp);

      // Set sizes
      this->mPadSize = setup.padSize();

      // Initialise temporary storage
      blockSize = this->mEBlockSize;
      Matrix tmpF = Matrix::Zero(fwdSize, blockSize);
      Matrix tmpB = Matrix::Zero(bwdSize, blockSize);

      // Create the spectral to physical plan
      const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
      this->mPlan = fftw_plan_many_r2r(1, fftSize, blockSize, tmpB.data(), NULL, 1, bwdSize, tmpF.data(), NULL, 1, fwdSize, bwdKind, Library::planFlag());
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
      this->mOddPlan = fftw_plan_many_r2r(1, fftSize, blockSize, tmpF.data(), NULL, 1, fwdSize, tmpB.data(), NULL, 1, bwdSize, ofwdKind, Library::planFlag());
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

   int WorlandProjector::lSize(const int l) const
   {
      return this->mWSize - l/2;
   }

   int WorlandProjector::lSize(const int i, const int ) const
   {
      return this->mWSize - i;
   }

   void WorlandProjector::io(const bool isEven) const
   {
      this->setPlan(isEven);

      this->io(this->mOutTmp.at(0).data(), this->mInTmp.at(0).data());
   }

   void WorlandProjector::input(const Matrix& in, const bool isEven, const bool needPadding) const
   {
      Matrix& inTmp = this->mInTmp.at(0);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven))
      {
         inTmp.block(0, start, this->mSpecSize, std::get<2>(loc)) = in.block(0, std::get<1>(loc), this->mSpecSize, std::get<2>(loc));
         start += std::get<2>(loc);
         std::get<3>(loc) = this->mSpecSize;
      }

      // Apply padding if required
      if(needPadding)
      {
         this->applyPadding(inTmp);
      }
   }

   void WorlandProjector::input(const MatrixZ& in, const bool isEven, const bool useReal, const bool needPadding) const
   {
      Matrix& inTmp = this->mInTmp.at(0);
      int start = 0;

      if(useReal)
      {
         for(auto& loc: *this->pLoc(isEven))
         {
            inTmp.block(0, start, this->mSpecSize, std::get<2>(loc)) = in.block(0, std::get<1>(loc), this->mSpecSize, std::get<2>(loc)).real();
            start += std::get<2>(loc);
            std::get<3>(loc) = this->mSpecSize;
         }
      } else
      {
         for(auto& loc: *this->pLoc(isEven))
         {
            inTmp.block(0, start, this->mSpecSize, std::get<2>(loc)) = in.block(0, std::get<1>(loc), this->mSpecSize, std::get<2>(loc)).imag();
            start += std::get<2>(loc);
            std::get<3>(loc) = this->mSpecSize;
         }
      }

      // Apply padding if required
      if(needPadding)
      {
         this->applyPadding(inTmp);
      }
   }

   void WorlandProjector::output(Matrix& rOut, const bool isEven) const
   {
      Matrix& outTmp = this->mOutTmp.at(0);
      int start = 0;
      for(auto loc: *this->pLoc(isEven))
      {
         rOut.block(0, std::get<1>(loc), rOut.rows(), std::get<2>(loc)) = outTmp.block(0, start, rOut.rows(), std::get<2>(loc));
         start += std::get<2>(loc);
      }

      // Zero unused values
      for(auto loc: this->mZLoc)
      {
         int s = std::get<1>(loc);
         int cols = std::get<2>(loc);
         rOut.block(0, s, rOut.rows(), cols).setZero();
      }
   }

   void WorlandProjector::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      Matrix& outTmp = this->mOutTmp.at(0);
      int start = 0;
      if(useReal)
      {
         for(auto loc: *this->pLoc(isEven))
         {
            rOut.block(0, std::get<1>(loc), rOut.rows(), std::get<2>(loc)).real() = outTmp.block(0, start, rOut.rows(), std::get<2>(loc));
            start += std::get<2>(loc);
         }

         // Zero unused values
         for(auto loc: this->mZLoc)
         {
            int s = std::get<1>(loc);
            int cols = std::get<2>(loc);
            rOut.block(0, s, rOut.rows(), cols).setZero();
         }
      } else
      {
         for(auto loc: *this->pLoc(isEven))
         {
            rOut.block(0, std::get<1>(loc), rOut.rows(), std::get<2>(loc)).imag() = outTmp.block(0, start, rOut.rows(), std::get<2>(loc));
            start += std::get<2>(loc);
         }

         // Zero unused values
         for(auto loc: this->mZLoc)
         {
            int s = std::get<1>(loc);
            int cols = std::get<2>(loc);
            rOut.block(0, s, rOut.rows(), cols).setZero();
         }
      }
   }

   void WorlandProjector::applyPadding(Matrix& rData, const int extraRows) const
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

   void WorlandProjector::applyFft() const
   {
      fftw_execute_r2r(*this->mpPlan, const_cast<MHDFloat *>(this->mpIn), this->mpOut);
   }

   void WorlandProjector::partialBackwardWorland(const int l, const int i0, const int start, const int cols, const std::vector<Matrix>& banded, Matrix& in) const
   {
      for(int i = i0; i >= l/2; --i)
      {
         this->applyPair(in, start, cols, this->lSize(i,l), banded.at(i));
      }
   }

   void WorlandProjector::backwardWorland(const bool isEven, const unsigned int id) const
   {
      Matrix& inTmp = this->mInTmp.at(id);
      const std::vector<Matrix>& banded = this->banded(isEven);

      std::vector<LocationType>* pL = this->pLoc(isEven,id);
      if(pL->begin() != pL->end())
      {
         int cols = 0;
         int maxL = std::get<4>(*pL->rbegin());
         int start = this->blockSize(isEven);
         int l;
         for(auto loc = pL->rbegin(); loc != pL->rend(); ++loc)
         {
            l = std::get<4>(*loc);
            this->partialBackwardWorland(l, maxL/2-1, start, cols, banded, inTmp);
            maxL = l;
            cols += std::get<2>(*loc);
            start -= std::get<2>(*loc);
         }

         // Last iteration if l=0/l=1 not in list
         if(std::get<4>(pL->at(0)) > 1)
         {
            l = static_cast<int>(!this->isPhysEven(isEven));
            this->partialBackwardWorland(l, maxL/2-1, start, cols, banded, inTmp);
         }

         // Rescale first mode for l = 0 for FFT
         if(std::get<4>(pL->at(0)) == 0)
         {
            inTmp.block(0, 0, 1, std::get<2>(pL->at(0))) *= std::sqrt(2.0);
         }

         // reset current l
         this->resetLocations(isEven, id);
      }
   }

   void WorlandProjector::lowerAlpha(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      Matrix& wTmp = this->workTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         Matrix U;
         this->buildShiftU(U, this->lSize(l), alpha, l-0.5, norm);
         this->applyTriSolve(wTmp, start, cols, U);
         start += cols;
      }
   }

   void WorlandProjector::lowerBeta(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      Matrix& wTmp = this->workTmp(id);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l > 0)
         {
            Matrix V;
            this->buildShiftV(V, std::get<3>(loc), alpha, l-0.5, norm);
            this->applyTriSolve(wTmp, start, cols, V);
            if(l == 1)
            {
               wTmp.block(0, start, 1, std::get<2>(loc)).array() *= 1.0/std::sqrt(2.0);
            }
            std::get<4>(loc)--;
         } else
         {
            wTmp.block(0, start, this->lSize(l), cols) /= norm;
            std::get<4>(loc)++;
         }
         start += cols;
      }
   }

   void WorlandProjector::lowerR2Beta(const MHDFloat alpha, const bool isEven, const unsigned int id, const MHDFloat norm) const
   {
      Matrix& wTmp = this->workTmp(id);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l < 0)
         {
            wTmp.block(0, start, wTmp.rows(), cols).setZero();
         } else if(l > 1)
         {
            Matrix M;
            this->buildShiftM(M, this->lSize(l), alpha, l-0.5, norm);
            this->applyTriProduct(wTmp, start, cols, M);
            std::get<4>(loc)--;
         } else
         {
            wTmp.block(0, start, this->lSize(l), cols) *= norm;
            std::get<4>(loc)--;
         }
         start += cols;
      }
   }

   void WorlandProjector::initBanded(const bool isEven) const
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
         int maxL = std::get<4>(*this->pLoc(isEven,0)->rbegin());
         int l = static_cast<int>(!this->isPhysEven(isEven));
         //for(int i = maxL/2-1; i >= l/2; --i)
         for(int i = l/2; i <= maxL/2-1; ++i)
         {
            Matrix PS;
            this->buildShiftPair(PS, false, i, this->lSize(i,l), J, normV, normM, false);
            pBanded->push_back(PS);
         }
      }
   }

}
}
}
}
}
