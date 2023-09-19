/**
 * @file WorlandProjector.cu
 * @brief Source of the interface for a generic cuFFT based Worland projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//
#include "Eigen/src/misc/blas.h"

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/WorlandProjector.hpp"

// Project includes
//
#include "Types/Constants.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuDctTools.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuMatrix.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   WorlandProjector::WorlandProjector()
      : IWorlandBackend(1,1)
   {
   }

   WorlandProjector::~WorlandProjector()
   {
   }

   void WorlandProjector::init(const SetupType& setup, const int lshift, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      //Initialize parent
      IWorlandBackend::init(setup, lshift, lshiftOnlyParity, alwaysZeroNegative);
      this->mpcuWorkTmp = &this->mcuInTmp;

      // Set sizes
      this->mPadSize = setup.padSize();

      // Set internal spectral resolution
      this->setWSize();

      // Initialize temporary storage
      int blockSize = std::max(this->mEBlockSize,this->mOBlockSize);
      // Input temporary storage
      this->initStorage(this->mBwdSize, blockSize, 2, this->mcuInTmp);
      // Output temporary storage
      this->initStorage(this->mFwdSize, blockSize, 1, this->mcuOutTmp);

      // Initialize Jacobi shift matrices
      this->initJ();

      // Create the complex to real plan
      this->mPlans.reserve(2*this->mNBatches);
      this->mcuFwd.reserve(this->mNBatches);
      this->mcuBwd.reserve(this->mNBatches);
      this->mcuWork.reserve(this->mNBatches);
      this->mcuWorkZ.reserve(this->mNBatches);

      // Initialize EVEN cuFFT plans and data
      int cuBwdSize = this->mFwdSize/2 + 1;
      int fwdSize = this->mFwdSize;
      int  *fftSize = &fwdSize;
      int cuBlockSize = blockSize/this->mNBatches + (blockSize%this->mNBatches > 0);
      int cuBatch = this->mEBlockSize/this->mNBatches +( this->mEBlockSize%this->mNBatches > 0);
      for(int i = 0; i < this->mNBatches; i++)
      {
         // Initialise GPU storage
         this->mcuFwd.push_back(0);
         this->mcuBwd.push_back(0);
         this->mcuWork.push_back(0);
         cudaMalloc((void**)&(this->mcuBwd.back()), sizeof(cufftDoubleComplex)*cuBwdSize*cuBlockSize);
         cudaMalloc((void**)&(this->mcuFwd.back()), sizeof(cufftDoubleReal)*this->mFwdSize*cuBlockSize);
         cudaMalloc((void**)&(this->mcuWork.back()), sizeof(cufftDoubleReal)*this->mFwdSize*cuBlockSize);

         // Create cuFFT plan
         cufftHandle plan;
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, cuBatch) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.at(i));
      }

      // Initialize ODD cuFFT plans and data
      fwdSize = this->mFwdSize/2;
      fftSize = &fwdSize;
      cuBatch = this->mOBlockSize/this->mNBatches +( this->mOBlockSize%this->mNBatches > 0);
      for(int i = 0; i < this->mNBatches; i++)
      {
         // Initialise GPU storage
         this->mcuWorkZ.push_back(0);
         cudaMalloc((void**)&(this->mcuWorkZ.back()), sizeof(cufftDoubleComplex)*fwdSize*cuBlockSize);

         // Create cuFFT plan
         cufftHandle plan;
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2Z, cuBatch) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.at(i));
      }
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

      this->io(this->mcuOutTmp.at(0).data(), this->mcuInTmp.at(0).data());
   }

   void WorlandProjector::input(const Matrix& in, const bool isEven, const bool needPadding) const
   {
      GpuMatrix& cuIn = this->mcuInTmp.at(0);
      int start = 0;
      for(auto& loc: *this->pLoc(isEven))
      {
         cublasSetMatrix(this->mSpecSize, std::get<2>(loc), sizeof(double), in.data()+std::get<1>(loc)*in.rows(), in.rows(), cuIn.data()+start*cuIn.rows(), cuIn.rows());
         start += std::get<2>(loc);
         std::get<3>(loc) = this->mSpecSize;
      }

      // Apply padding if required
      if(needPadding)
      {
         this->applyPadding(cuIn);
      }
   }

   void WorlandProjector::input(const MatrixZ& in, const bool isEven, const bool useReal, const bool needPadding) const
   {
      GpuMatrix& cuIn = this->mcuInTmp.at(0);
      int start = 0;

      int zShift;
      if(useReal)
      {
         zShift = 0;
      } else
      {
         zShift = 1;
      }

      for(auto& loc: *this->pLoc(isEven))
      {
         for(int k = 0; k < std::get<2>(loc); k++)
         {
            cublasSetVector(this->mSpecSize, sizeof(double), reinterpret_cast<const double*>(in.data())+2*(std::get<1>(loc)+k)*in.rows()+zShift, 2, cuIn.data()+(start+k)*cuIn.rows(), 1);
         }
         start += std::get<2>(loc);
         std::get<3>(loc) = this->mSpecSize;
      }

      // Apply padding if required
      if(needPadding)
      {
         this->applyPadding(cuIn);
      }
   }

   void WorlandProjector::output(Matrix& rOut, const bool isEven) const
   {
      GpuMatrix& cuOut = this->mcuOutTmp.at(0);
      int start = 0;
      for(auto loc: *this->pLoc(isEven))
      {
         for(int k = 0; k < std::get<2>(loc); k++)
         {
            cublasGetVector(rOut.rows(), sizeof(double), cuOut.data()+(start+k)*cuOut.rows(), 1, rOut.data()+(std::get<1>(loc)+k)*rOut.rows(), 1);
         }
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
      GpuMatrix& cuOut = this->mcuOutTmp.at(0);
      int start = 0;
      if(useReal)
      {
         for(auto loc: *this->pLoc(isEven))
         {
            int zShift = 0;
            for(int k = 0; k < std::get<2>(loc); k++)
            {
               cublasGetVector(rOut.rows(), sizeof(double), cuOut.data()+(start+k)*cuOut.rows(), 1, reinterpret_cast<double*>(rOut.data())+2*(std::get<1>(loc)+k)*rOut.rows()+zShift, 2);
            }
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
         int zShift = 1;
         for(auto loc: *this->pLoc(isEven))
         {
            for(int k = 0; k < std::get<2>(loc); k++)
            {
               cublasGetVector(rOut.rows(), sizeof(double), cuOut.data()+(start+k)*cuOut.rows(), 1, reinterpret_cast<double*>(rOut.data())+2*(std::get<1>(loc)+k)*rOut.rows()+zShift, 2);
            }
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

   void WorlandProjector::applyPadding(GpuMatrix& rData, const int extraRows) const
   {
      double zero = 0.0;
      int zRows;
      if(extraRows >= 0)
      {
         zRows = this->mPadSize-extraRows;
      } else
      {
         zRows = -extraRows;
      }
      for(int k = 0; k < zRows; k++)
      {
         cublasDscal(this->mHBlas, rData.cols(), &zero, rData.data() + rData.rows()-zRows+k, rData.rows());
      }
   }

   void WorlandProjector::applyEvenFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      dim3 kblocks(16,16);
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNBatches;
         int blocks = this->mcuEBlockSize.at(i);
         Gpu::buildIDCTInput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuBwd.at(sid),this->mpIn+bshift, this->mBwdSize, blocks);
         bshift += this->mBwdSize*blocks;
         cufftExecZ2D(this->mPlans.at(sid), this->mcuBwd.at(sid), this->mcuWork.at(sid));
         Gpu::buildIDCTOutput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mpOut+fshift,this->mcuWork.at(sid), this->mFwdSize, blocks);
         fshift += this->mFwdSize*blocks;
      }

      for(int i = 0; i < this->mNBatches; i++)
      {
         if(cudaStreamSynchronize(this->mStreams.at(i)) != cudaSuccess)
         {
            throw std::logic_error("CUFFT Error: Synchronization failed");
         }
      }
   }

   void WorlandProjector::applyOddFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      dim3 kblocks(16,16);
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNBatches;
         int blocks = this->mcuOBlockSize.at(i);
         Gpu::buildDCT4Input<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuBwd.at(sid),this->mpIn + bshift, this->mBwdSize, blocks);
         bshift += this->mBwdSize*blocks;
         cufftExecZ2Z(this->mPlans.at(sid+this->mNBatches), this->mcuBwd.at(sid), this->mcuWorkZ.at(sid), CUFFT_FORWARD);
         Gpu::buildDCT4Output<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mpOut+fshift,this->mcuWorkZ.at(sid), this->mFwdSize, blocks, 2.0);
         fshift += this->mFwdSize*blocks;
      }

      for(int i = 0; i < this->mNBatches; i++)
      {
         if(cudaStreamSynchronize(this->mStreams.at(i)) != cudaSuccess)
         {
            throw std::logic_error("CUFFT Error: Synchronization failed");
         }
      }
   }

   void WorlandProjector::applyFft() const
   {
     if(this->mPlanIsEven)
     {
       this->applyEvenFft();
     } else
     {
       this->applyOddFft();
     }
   }

   void WorlandProjector::partialBackwardWorland(const int l, const int i0, const int start, const int cols, const GpuMatrix& J, double* in, const int iSize, const double normV, const double normM) const
   {
      int sId = 0;
      for(int i = i0; i >= l/2; --i)
      {
         this->mcuV.reshape(GpuMatrix::KEEP_SIZE, this->lSize(i,l));
         this->mcuM.reshape(GpuMatrix::KEEP_SIZE, this->lSize(i,l));
         this->buildShiftPair(sId, this->mcuV, this->mcuM, i, J, normV, normM, false);
         this->applyPair(sId, in, iSize, start, cols, this->mcuM, this->mcuV);
         this->synchronize(sId);
      }
   }

   void WorlandProjector::backwardWorland(const bool isEven, const unsigned int id) const
   {
      const GpuMatrix& J = this->J(isEven);
      GpuMatrix& cuIn = this->cuWorkTmp(id);

      std::vector<LocationType>* pL = this->pLoc(isEven,id);
      int cols = 0;
      int maxL = std::get<4>(*pL->rbegin());
      int start = this->blockSize(isEven);
      int l;
      const double normV = 1.0;
      const double normM = 1.0;
      for(auto loc = pL->rbegin(); loc != pL->rend(); ++loc)
      {
         l = std::get<4>(*loc);
         this->partialBackwardWorland(l, maxL/2-1, start, cols, J, cuIn.data(), cuIn.rows(), normV, normM);
         maxL = l;
         cols += std::get<2>(*loc);
         start -= std::get<2>(*loc);
      }

      // Last iteration if l=0/l=1 not in list
      if(std::get<4>(pL->at(0)) > 1)
      {
         l = static_cast<int>(!this->isPhysEven(isEven));
         this->partialBackwardWorland(l, maxL/2-1, start, cols, J, cuIn.data(), cuIn.rows(), normV, normM);
      }

      // Rescale first mode for l = 0 for FFT
      if(std::get<4>(pL->at(0)) == 0)
      {
         double scale;
         scale = std::sqrt(2.0);
         cublasDscal(this->mHBlas, std::get<2>(pL->at(0)), &scale, cuIn.data(), cuIn.rows());
      }

      // reset current l
      this->resetLocations(isEven, id);
   }

   void WorlandProjector::lowerAlpha(const double alpha, const bool isEven, const unsigned int id, const double norm) const
   {
      GpuMatrix& cuIn = this->cuWorkTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         this->mcuU.reshape(GpuMatrix::KEEP_SIZE, this->lSize(l));
         this->buildShiftU(this->mcuU, alpha, l-0.5, norm);
         this->applyTriSolve(cuIn.data(), cuIn.rows(), start, cols, this->mcuU);
         start += cols;
      }
   }

   void WorlandProjector::lowerBeta(const double alpha, const bool isEven, const unsigned int id, const double norm) const
   {
      GpuMatrix& cuIn = this->cuWorkTmp(id);

      int start = 0;
      double scale;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l > 0)
         {
            this->mcuV.reshape(GpuMatrix::KEEP_SIZE, std::get<3>(loc));
            this->buildShiftV(this->mcuV, alpha, l-0.5, norm);
            this->applyTriSolve(cuIn.data(), cuIn.rows(), start, cols, this->mcuV);
            if(l == 1)
            {
               scale = 1.0/std::sqrt(2.0);
               cublasDscal(this->mHBlas, std::get<2>(loc), &scale, cuIn.data()+start*cuIn.rows(), cuIn.rows());
            }
            std::get<4>(loc)--;
         } else
         {
            scale = 1.0/norm;
            for(int k = 0; k < cols; k++)
            {
               cublasDscal(this->mHBlas, this->lSize(l), &scale, cuIn.data()+(start+k)*cuIn.rows(), 1);
            }
            std::get<4>(loc)++;
         }
         start += cols;
      }
   }

   void WorlandProjector::lowerR2Beta(const double alpha, const bool isEven, const unsigned int id, const double norm) const
   {
      GpuMatrix& cuIn = this->cuWorkTmp(id);
      int start = 0;
      double scale;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l < 0)
         {
            scale = 0.0;
            cublasDscal(this->mHBlas, cuIn.rows()*cols, &scale, cuIn.data()+start*cuIn.rows(), 1);
         } else if(l > 1)
         {
            this->mcuM.reshape(GpuMatrix::KEEP_SIZE, this->lSize(l));
            this->buildShiftM(this->mcuM, alpha, l-0.5, norm);
            this->applyTriProduct(cuIn.data(), cuIn.rows(), start, cols, this->mcuM);
            std::get<4>(loc)--;
         } else
         {
            scale = norm;
            for(int k = 0; k < cols; k++)
            {
               cublasDscal(this->mHBlas, this->lSize(l), &scale, cuIn.data()+(start+k)*cuIn.rows(), 1);
            }
            std::get<4>(loc)--;
         }
         start += cols;
      }
   }

}
}
}
}
}
