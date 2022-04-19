/**
 * @file ChebyshevProjector.cu
 * @brief Source of the interface for a generic cuFFT based Chebyshev projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/ChebyshevProjector.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Backend/CuFft/GpuDctTools.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   ChebyshevProjector::ChebyshevProjector()
      : IChebyshevBackend(1,1)
   {
   }

   ChebyshevProjector::~ChebyshevProjector()
   {
   }

   void ChebyshevProjector::init(const SetupType& setup) const
   {
      //Initialize parent
      IChebyshevBackend::init(setup);

      // Set sizes
      int fwdSize = this->mFwdSize;
      int cuBwdSize = this->mFwdSize/2 + 1;
      int  *fftSize = &fwdSize;
      this->mPadSize = setup.padSize();

      // Initialize storage
      this->mTmp.setZero(this->mBwdSize, this->mBlockSize);
      this->mTmpComp.setZero(this->mBwdSize, this->mBlockSize);

      // Create the complex to real plan
      this->mPlans.reserve(this->mNStreams);
      this->mStreams.reserve(this->mNStreams);
      this->mcuFwd.reserve(this->mNStreams);
      this->mcuBwd.reserve(this->mNStreams);
      this->mcuWork.reserve(this->mNStreams);

      for(int i = 0; i < this->mNStreams; i++)
      {
         // Initialise GPU storage
         this->mcuFwd.push_back(0);
         this->mcuBwd.push_back(0);
         this->mcuWork.push_back(0);
         cudaMalloc((void**)&(this->mcuBwd.back()), sizeof(cufftDoubleComplex)*cuBwdSize*this->mBlockSize);
         cudaMalloc((void**)&(this->mcuFwd.back()), sizeof(cufftDoubleReal)*this->mFwdSize*this->mBlockSize);
         cudaMalloc((void**)&(this->mcuWork.back()), sizeof(cufftDoubleReal)*this->mFwdSize*this->mBlockSize);

         // Create CUDA stream
         cudaStream_t stream;
         cudaStreamCreate(&stream);
         this->mStreams.push_back(stream);

         // Create cuFFT plan
         cufftHandle plan;
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, this->mBlockSize) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.back());
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
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      dim3 kblocks(16,16);
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNStreams;
         cudaMemcpyAsync(this->mcuWork.at(sid), this->mpIn+bshift, sizeof(cufftDoubleReal)*this->mBwdSize*this->mBlockSize,cudaMemcpyHostToDevice, this->mStreams.at(sid));
         bshift += this->mBwdSize*this->mBlockSize;
         Gpu::buildIDCTInput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuBwd.at(sid),this->mcuWork.at(sid), this->mBwdSize, this->mBlockSize);
         cufftExecZ2D(this->mPlans.at(sid), this->mcuBwd.at(sid), this->mcuWork.at(sid));
         Gpu::buildIDCTOutput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuFwd.at(sid),this->mcuWork.at(sid), this->mFwdSize, this->mBlockSize);
         cudaMemcpyAsync(this->mpOut+fshift, this->mcuFwd.at(sid), sizeof(cufftDoubleReal)*this->mFwdSize*this->mBlockSize,cudaMemcpyDeviceToHost, this->mStreams.at(sid));
         fshift += this->mFwdSize*this->mBlockSize;
      }

      for(int i = 0; i < this->mNStreams; i++)
      {
         if(cudaStreamSynchronize(this->mStreams.at(i)) != cudaSuccess)
         {
            throw std::logic_error("CUFFT Error: Synchronization failed");
         }
      }
   }

   void ChebyshevProjector::addSolver(const int extraRows) const
   {
      this->mspSolver = std::make_shared<Fftw::DifferentialSolver>(this->mSpecSize, this->mBlockSize, extraRows);
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

   Fftw::DifferentialSolver& ChebyshevProjector::solver() const
   {
      return *this->mspSolver;
   }

}
}
}
}
}
