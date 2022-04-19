/**
 * @file MixedProjector.cu
 * @brief Source of the interface for a generic cuFFT based mixed projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/MixedProjector.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   MixedProjector::MixedProjector()
      : IMixedBackend(5, 10)
   {
   }

   MixedProjector::~MixedProjector()
   {
   }

   void MixedProjector::init(const SetupType& setup) const
   {
      //Initialize parent
      IMixedBackend::init(setup);

      int fwdSize = this->mFwdSize;
      int  *fftSize = &fwdSize;
      this->mPadSize = setup.padSize();

      // Initialise temporary storage
      this->mTmp.setZero(this->mBwdSize, setup.blockSize());

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
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2D, this->mBlockSize) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.back());
      }
   }

   void MixedProjector::input(const MatrixZ& in) const
   {
      this->mTmp.topRows(this->mSpecSize) = in.topRows(this->mSpecSize);

      this->applyPadding(this->mTmp);
   }

   void MixedProjector::inputDiff(const MatrixZ& in, const int order, const double scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*Math::cI;
         this->mTmp.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mSpecSize);
      } else
      {
         double sgn = std::pow(-1.0,(order/2)%2)*scale;
         this->mTmp.topRows(this->mSpecSize) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mSpecSize);
      }

      this->applyPadding(this->mTmp);
   }

   void MixedProjector::output(double* out) const
   {
      this->io(out, this->mTmp.data());
   }

   void MixedProjector::io(double* out, const MHDComplex* in) const
   {
      this->mpOut = out;
      this->mpIn = in;
   }

   void MixedProjector::applyPadding(MatrixZ& rData) const
   {
      // Set the m=0 values to zero
      rData.row(0).imag().setConstant(0);

      // Set the padded values to zero
      rData.bottomRows(this->mPadSize).setZero();
   }

   void MixedProjector::applyFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNStreams;
         cudaMemcpyAsync(this->mcuBwd.at(sid), this->mpIn+bshift, sizeof(cufftDoubleComplex)*this->mBwdSize*this->mBlockSize,cudaMemcpyHostToDevice, this->mStreams.at(sid));
         bshift += this->mBwdSize*this->mBlockSize;
         cufftExecZ2D(this->mPlans.at(sid), this->mcuBwd.at(sid), this->mcuFwd.at(sid));
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

}
}
}
}
}
