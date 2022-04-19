/**
 * @file ChebyshevIntegrator.cu
 * @brief Source of the interface for a generic cuFFT based Chebyshev integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/ChebyshevIntegrator.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuDctTools.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   ChebyshevIntegrator::ChebyshevIntegrator()
      : IChebyshevBackend(1,1), mOutMap(NULL,0,0)
   {
   }

   ChebyshevIntegrator::~ChebyshevIntegrator()
   {
   }

   void ChebyshevIntegrator::init(const SetupType& setup) const
   {
      //Initialize parent
      IChebyshevBackend::init(setup);

      int fwdSize = this->mFwdSize;
      int cuBwdSize = this->mFwdSize/2 + 1;
      int  *fftSize = &fwdSize;

      // Set transform scaling
      this->mFftScaling = 1.0/static_cast<double>(this->mFwdSize);

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
         // Create the complex to real plan
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, this->mBlockSize) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.back());
      }

      // Initialize temporary storage
      this->mTmp.setZero(this->mFwdSize, this->mBlockSize);
   }

   void ChebyshevIntegrator::io(double* out, const double* in) const
   {
      IChebyshevBackend::io(out, in);

      new (&this->mOutMap) Eigen::Map<Matrix>(this->mpOut, this->mBwdSize, this->mBlockSize);
   }

   void ChebyshevIntegrator::io(Matrix& rOut, const Matrix& in) const
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

   void ChebyshevIntegrator::input(const MatrixZ& in, const bool useReal) const
   {
      if(useReal)
      {
         this->mTmpComp = in.real();
      } else
      {
         this->mTmpComp = in.imag();
      }

      this->io(this->mTmp.data(), this->mTmpComp.data());
   }

   void ChebyshevIntegrator::applyFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      dim3 kblocks(16,16);
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNStreams;
         cudaMemcpyAsync(this->mcuFwd.at(sid), this->mpIn+fshift, sizeof(cufftDoubleReal)*this->mFwdSize*this->mBlockSize,cudaMemcpyHostToDevice, this->mStreams.at(sid));
         fshift += this->mFwdSize*this->mBlockSize;
         Gpu::buildDCTInput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuWork.at(sid),this->mcuFwd.at(sid), this->mFwdSize, this->mBlockSize);
         cufftExecD2Z(this->mPlans.at(sid), this->mcuWork.at(sid), this->mcuBwd.at(sid));
         Gpu::buildDCTOutput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuWork.at(sid),this->mcuBwd.at(sid), this->mBwdSize, this->mBlockSize);
         cudaMemcpyAsync(this->mpOut+bshift, this->mcuWork.at(sid), sizeof(cufftDoubleReal)*this->mBwdSize*this->mBlockSize,cudaMemcpyDeviceToHost, this->mStreams.at(sid));
         bshift += this->mBwdSize*this->mBlockSize;
      }

      for(int i = 0; i < this->mNStreams; i++)
      {
         if(cudaStreamSynchronize(this->mStreams.at(i)) != cudaSuccess)
         {
            throw std::logic_error("CUFFT Error: Synchronization failed");
         }
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
      rOut.topRows(this->mSpecSize) = this->mFftScaling*this->mOutMap.topRows(this->mSpecSize);
   }

   void ChebyshevIntegrator::outputSpectral(Matrix& rOut) const
   {
      if(this->mMeanOp.size() > 0)
      {
         rOut.block(0, 0, this->mSpecSize, 1) = this->mFftScaling*this->mMeanOp*this->mOutMap.block(0,0,this->mMeanOp.cols(),1);
         rOut.block(0, 1, this->mSpecSize, rOut.cols()-1) = this->mFftScaling*this->mSpecOp*this->mOutMap.block(0,1,this->mSpecOp.cols(), rOut.cols()-1);
      } else
      {
         rOut.topRows(this->mSpecSize) = this->mFftScaling*this->mSpecOp*this->mOutMap.topRows(this->mSpecOp.cols());
      }
   }

   void ChebyshevIntegrator::output(MatrixZ& rOut, const bool useReal) const
   {
      if(this->mMeanOp.size() > 0)
      {
         if(useReal)
         {
            rOut.block(0, 0, this->mSpecSize, 1).real() = this->mFftScaling*this->mMeanOp*this->mOutMap.block(0,0,this->mMeanOp.cols(),1);
            rOut.block(0, 1, this->mSpecSize, rOut.cols()-1).real() = this->mFftScaling*this->mSpecOp*this->mOutMap.block(0,1,this->mSpecOp.cols(), rOut.cols()-1);
         } else
         {
            rOut.block(0, 0, this->mSpecSize, 1).imag() = this->mFftScaling*this->mMeanOp*this->mOutMap.block(0,0,this->mMeanOp.cols(),1);
            rOut.block(0, 1, this->mSpecSize, rOut.cols()-1).imag() = this->mFftScaling*this->mSpecOp*this->mOutMap.block(0,1,this->mSpecOp.cols(), rOut.cols()-1);
         }
      } else
      {
         if(useReal)
         {
            rOut.topRows(this->mSpecSize).real() = this->mFftScaling*this->mOutMap.topRows(this->mSpecSize);
         } else
         {
            rOut.topRows(this->mSpecSize).imag() = this->mFftScaling*this->mOutMap.topRows(this->mSpecSize);
         }
      }
   }

   void ChebyshevIntegrator::outputSpectral(MatrixZ& rOut, const bool useReal) const
   {
      if(useReal)
      {
         rOut.topRows(this->mSpecSize).real() = this->mFftScaling*this->mSpecOp*this->mOutMap.topRows(this->mSpecOp.cols());
      } else
      {
         rOut.topRows(this->mSpecSize).imag() = this->mFftScaling*this->mSpecOp*this->mOutMap.topRows(this->mSpecOp.cols());
      }
   }

}
}
}
}
}
