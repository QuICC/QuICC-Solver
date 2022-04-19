/**
 * @file ChebyshevEnergy.cu
 * @brief Source of the interface for a generic cuFFT based Chebyshev energy reductor
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/ChebyshevEnergy.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuDctTools.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   ChebyshevEnergy::ChebyshevEnergy()
      : IChebyshevBackend(1,1)
   {
   }

   ChebyshevEnergy::~ChebyshevEnergy()
   {
   }

   void ChebyshevEnergy::init(const SetupType& setup) const
   {
      //Initialize parent
      IChebyshevBackend::init(setup);

      this->mNrgFwdSize = 2*setup.fwdSize()+4;
      this->mNrgBwdSize = 2*setup.bwdSize()+4;
      int cuBwdSize = this->mNrgFwdSize/2 + 1;
      int fwdSize = this->mNrgFwdSize;
      int  *fftSize = &fwdSize;

      // Set transform scaling
      this->mFftScaling = 1.0/static_cast<double>(this->mNrgFwdSize);

      // Initialize storage
      this->mTmp.setZero(this->mNrgBwdSize, this->mBlockSize);
      this->mTmpComp.setZero(this->mNrgBwdSize, this->mBlockSize);
      this->mTmpMid.setZero(this->mNrgFwdSize, this->mBlockSize);

      // Compute energy weights
      this->computeEWeights(this->mNrgBwdSize, setup.lower(), setup.upper());

      // Initialise temporary storage
      Matrix tmpF = Matrix::Zero(this->mNrgFwdSize, this->mBlockSize);
      Matrix tmpB = Matrix::Zero(this->mNrgBwdSize, this->mBlockSize);

      // Create the complex to real plan
      this->mPlans.reserve(2*this->mNStreams);
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
         cudaMalloc((void**)&(this->mcuFwd.back()), sizeof(cufftDoubleReal)*this->mNrgFwdSize*this->mBlockSize);
         cudaMalloc((void**)&(this->mcuWork.back()), sizeof(cufftDoubleReal)*this->mNrgFwdSize*this->mBlockSize);

         // Create CUDA stream
         cudaStream_t stream;
         cudaStreamCreate(&stream);
         this->mStreams.push_back(stream);

         // Create the complex to rea; cuFFT plan
         cufftHandle plan;
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, this->mBlockSize) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.back());
      }

      for(int i = 0; i < this->mNStreams; i++)
      {

         // Create the real to complex plan
         cufftHandle plan;
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, this->mBlockSize) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.at(i));
      }
   }

   void ChebyshevEnergy::computeEWeights(const int size, const double lower, const double upper) const
   {
      // Initialize energy weights
      this->mEWeights = Array::Zero(size);
      double a = (upper - lower)/2.0;
      for(int i = 0; i < size/2; i++)
      {
         double n = 2.0*i;
         this->mEWeights(2*i) = 2.0*a*(2.0/(1.0 - n*n));
      }
      this->mEWeights(0) *= 0.5;
   }

   void ChebyshevEnergy::applyPadding(Matrix& rData, const int extraRows) const
   {
      // Set the padded values to zero
      rData.bottomRows(this->mFwdSize+this->mPadSize-extraRows).setZero();
   }

   void ChebyshevEnergy::applyFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      dim3 kblocks(16,16);
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNStreams;
         cudaMemcpyAsync(this->mcuWork.at(sid), this->mpIn+bshift, sizeof(cufftDoubleReal)*this->mNrgBwdSize*this->mBlockSize,cudaMemcpyHostToDevice, this->mStreams.at(sid));
         bshift += this->mNrgBwdSize*this->mBlockSize;
         Gpu::buildIDCTInput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuBwd.at(sid),this->mcuWork.at(sid), this->mNrgBwdSize, this->mBlockSize);
         cufftExecZ2D(this->mPlans.at(sid), this->mcuBwd.at(sid), this->mcuWork.at(sid));
         Gpu::buildIDCTOutput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuFwd.at(sid),this->mcuWork.at(sid), this->mNrgFwdSize, this->mBlockSize);
         cudaMemcpyAsync(this->mpOut+fshift, this->mcuFwd.at(sid), sizeof(cufftDoubleReal)*this->mNrgFwdSize*this->mBlockSize,cudaMemcpyDeviceToHost, this->mStreams.at(sid));
         fshift += this->mNrgFwdSize*this->mBlockSize;
      }

      for(int i = 0; i < this->mNStreams; i++)
      {
         if(cudaStreamSynchronize(this->mStreams.at(i)) != cudaSuccess)
         {
            throw std::logic_error("CUFFT Error: Synchronization failed");
         }
      }
   }

   void ChebyshevEnergy::applyFwdFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      dim3 kblocks(16,16);
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = (i % this->mNStreams);
         cudaMemcpyAsync(this->mcuFwd.at(sid), this->mpIn+fshift, sizeof(cufftDoubleReal)*this->mNrgFwdSize*this->mBlockSize,cudaMemcpyHostToDevice, this->mStreams.at(sid));
         fshift += this->mNrgFwdSize*this->mBlockSize;
         Gpu::buildDCTInput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuWork.at(sid),this->mcuFwd.at(sid), this->mNrgFwdSize, this->mBlockSize);
         cufftExecD2Z(this->mPlans.at(sid+this->mNStreams), this->mcuWork.at(sid), this->mcuBwd.at(sid));
         Gpu::buildDCTOutput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuWork.at(sid),this->mcuBwd.at(sid), this->mNrgBwdSize, this->mBlockSize);
         cudaMemcpyAsync(this->mpOut+bshift, this->mcuWork.at(sid), sizeof(cufftDoubleReal)*this->mNrgBwdSize*this->mBlockSize,cudaMemcpyDeviceToHost, this->mStreams.at(sid));
         bshift += this->mNrgBwdSize*this->mBlockSize;
      }

      for(int i = 0; i < this->mNStreams; i++)
      {
         if(cudaStreamSynchronize(this->mStreams.at(i)) != cudaSuccess)
         {
            throw std::logic_error("CUFFT Error: Synchronization failed");
         }
      }
   }

   void ChebyshevEnergy::input(const Matrix& in, const bool needPadding) const
   {
      this->mTmp.topRows(this->mSpecSize) = in.topRows(this->mSpecSize);

      // Apply padding if required
      if(needPadding)
      {
         this->applyPadding(this->mTmp);
      }
   }

   void ChebyshevEnergy::input(const MatrixZ& in, const bool useReal, const bool needPadding) const
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

   void ChebyshevEnergy::io() const
   {
      this->io(this->mTmpComp.data(), this->mTmp.data());
   }

   void ChebyshevEnergy::setSpectralOperator(const SparseMatrix& mat) const
   {
      this->mSpecOp = mat;
   }

   void ChebyshevEnergy::outputSpectral(Matrix& rOut) const
   {
      assert(this->mSpecOp.cols() <= this->mTmp.rows());
      assert(this->mSpecOp.rows() <= this->mEWeights.rows());
      rOut.transpose() = this->mFftScaling*this->mEWeights.topRows(this->mSpecOp.rows()).transpose()*this->mSpecOp*this->mTmp.topRows(2*this->mSpecSize);
   }

   void ChebyshevEnergy::output(Matrix& rOut) const
   {
      int rows = 2*this->mSpecSize;
      rOut.transpose() = this->mFftScaling*this->mEWeights.topRows(rows).transpose()*this->mTmp.topRows(rows);
   }

   void ChebyshevEnergy::square(const bool isFirst) const
   {
      if(isFirst)
      {
         this->mTmpMid = this->mTmpComp.array().pow(2);
      } else
      {
         this->mTmpMid.array() += this->mTmpComp.array().pow(2);
      }
      this->io(this->mTmp.data(), this->mTmpMid.data());
   }

   void ChebyshevEnergy::addSolver(const int extraRows) const
   {
      this->mspSolver = std::make_shared<Fftw::DifferentialSolver>(this->mSpecSize, this->mBlockSize, extraRows);
   }

   void ChebyshevEnergy::getSolution(const int zeroRows, const int extraRows) const
   {
      this->solver().solve(zeroRows, this->mTmp);
      this->applyPadding(this->mTmp, extraRows);
   }

   Fftw::DifferentialSolver& ChebyshevEnergy::solver() const
   {
      return *this->mspSolver;
   }

}
}
}
}
}
