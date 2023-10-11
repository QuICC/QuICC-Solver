/**
 * @file MixedIntegrator.cu
 * @brief Source of the interface for a generic cuFFT based mixed integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/MixedIntegrator.hpp"

// Project includes
//
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   MixedIntegrator::MixedIntegrator()
      : IMixedBackend(5, 10), mOutMap(NULL,0,0)
   {
   }

   MixedIntegrator::~MixedIntegrator()
   {
   }

   void MixedIntegrator::init(const SetupType& setup) const
   {
      //Initialize parent
      IMixedBackend::init(setup);

      int fwdSize = this->mFwdSize;
      int  *fftSize = &fwdSize;

      // Set transform scaling
      this->mFftScaling = 1.0/static_cast<double>(this->mFwdSize);

      // Create the complex to real plan
      this->mPlans.reserve(this->mNStreams);
      this->mStreams.reserve(this->mNStreams);
      this->mcuFwd.reserve(this->mNStreams);
      this->mcuBwd.reserve(this->mNStreams);

      for(int i = 0; i < this->mNStreams; i++)
      {
         // Initialise GPU storage
         this->mcuFwd.push_back(0);
         this->mcuBwd.push_back(0);
         cudaMalloc((void**)&(this->mcuBwd.back()), sizeof(cufftDoubleComplex)*this->mBwdSize*this->mBlockSize);
         cudaMalloc((void**)&(this->mcuFwd.back()), sizeof(cufftDoubleReal)*this->mFwdSize*this->mBlockSize);

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
   }

   void MixedIntegrator::output(MatrixZ& rOut) const
   {
      rOut.topRows(this->mSpecSize) = this->mFftScaling*this->mOutMap.topRows(this->mSpecSize);
   }

   void MixedIntegrator::outputDiff(MatrixZ& rOut, const int order, const double scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*this->mFftScaling*Math::cI;
         rOut.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*this->mOutMap.topRows(this->mSpecSize);
      } else
      {
         double sgn = std::pow(-1.0,(order/2)%2)*this->mFftScaling*scale;
         rOut.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*this->mOutMap.topRows(this->mSpecSize);
      }
   }

   void MixedIntegrator::outputDiff(MatrixZ& rOut, const int order, const double scale, const std::map<int,MHDComplex>& mod) const
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
         double sgn = std::pow(-1.0,(order/2)%2)*this->mFftScaling;
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
         this->mTmp.setConstant(-42424242);
         this->io(this->mTmp.data(), in.data());
      }
   }

   void MixedIntegrator::io(MHDComplex* out, const double* in) const
   {
      this->mpOut = out;
      this->mpIn = in;

      new (&this->mOutMap) Eigen::Map<MatrixZ>(this->mpOut, this->mBwdSize, this->mBlockSize);
   }

   void MixedIntegrator::applyFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNStreams;
         cudaMemcpyAsync(this->mcuFwd.at(sid), this->mpIn+fshift, sizeof(cufftDoubleReal)*this->mFwdSize*this->mBlockSize,cudaMemcpyHostToDevice, this->mStreams.at(sid));
         fshift += this->mFwdSize*this->mBlockSize;
         cufftExecD2Z(this->mPlans.at(sid), this->mcuFwd.at(sid), this->mcuBwd.at(sid));
         cudaMemcpyAsync(this->mpOut+bshift, this->mcuBwd.at(sid), sizeof(cufftDoubleComplex)*this->mBwdSize*this->mBlockSize,cudaMemcpyDeviceToHost, this->mStreams.at(sid));
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

}
}
}
}
}
