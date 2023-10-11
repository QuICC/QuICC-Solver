/**
 * @file ComplexProjector.cu
 * @brief Source of the interface for a generic cuFFt based complex projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/ComplexProjector.hpp"

// Project includes
//
#include "Types/Math.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   ComplexProjector::ComplexProjector()
      : IComplexBackend(5, 10)
   {
   }

   ComplexProjector::~ComplexProjector()
   {
   }

   void ComplexProjector::init(const SetupType& setup) const
   {
      //Initialize parent
      IComplexBackend::init(setup);

      int fwdSize = this->mFwdSize;
      int  *fftSize = &fwdSize;
      this->mPadSize = setup.padSize();

      // Initialise temporary storage
      this->mTmp.setZero(this->mBwdSize, this->mBlockSize);

      this->initMeanBlocks(setup.idBlocks());

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
         cudaMalloc((void**)&(this->mcuFwd.back()), sizeof(cufftDoubleComplex)*this->mFwdSize*this->mBlockSize);

         // Create CUDA stream
         cudaStream_t stream;
         cudaStreamCreate(&stream);
         this->mStreams.push_back(stream);

         // Create cuFFT plan
         cufftHandle plan;
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2Z, this->mBlockSize) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.back());
      }
   }

   void ComplexProjector::applyFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNStreams;
         cudaMemcpyAsync(this->mcuBwd.at(sid), this->mpIn+bshift, sizeof(cufftDoubleComplex)*this->mBwdSize*this->mBlockSize,cudaMemcpyHostToDevice, this->mStreams.at(sid));
         bshift += this->mBwdSize*this->mBlockSize;
         cufftExecZ2Z(this->mPlans.at(sid), this->mcuBwd.at(sid), this->mcuFwd.at(sid), CUFFT_INVERSE);
         cudaMemcpyAsync(this->mpOut+fshift, this->mcuFwd.at(sid), sizeof(cufftDoubleComplex)*this->mFwdSize*this->mBlockSize,cudaMemcpyDeviceToHost, this->mStreams.at(sid));
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

   void ComplexProjector::input(const MatrixZ& in) const
   {
      this->mTmp.topRows(this->mPosN) = in.topRows(this->mPosN);
      this->mTmp.bottomRows(this->mNegN) = in.bottomRows(this->mNegN);

      this->forceConjugate(this->mTmp);
      this->applyPadding(this->mTmp);
   }

   void ComplexProjector::inputMean(const MatrixZ& in) const
   {
      // Initialize to zero
      this->mTmp.setZero();

      // Set the mean
      for(auto it = this->mMeanBlocks.cbegin(); it != this->mMeanBlocks.cend(); ++it)
      {
         this->mTmp.block(0, it->first, 1, it->second) = in.block(0, it->first, 1, it->second);
      }
   }

   void ComplexProjector::inputDiff(const MatrixZ& in, const int order, const double scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*Math::cI;
         this->mTmp.topRows(this->mPosN) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mPosN);
         this->mTmp.bottomRows(this->mNegN) = (sgn*(scale*this->negativeK()).array().pow(order).matrix()).asDiagonal()*in.bottomRows(this->mNegN);
      } else
      {
         double sgn = std::pow(-1.0,(order/2)%2);
         this->mTmp.topRows(this->mPosN) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*in.topRows(this->mPosN);
         this->mTmp.bottomRows(this->mNegN) = (sgn*(scale*this->negativeK()).array().pow(order).matrix()).asDiagonal()*in.bottomRows(this->mNegN);
      }

      this->forceConjugate(this->mTmp);
      this->applyPadding(this->mTmp);
   }

   void ComplexProjector::inputDiff2D(const MatrixZ& in, const std::vector<std::pair<int,int> >& orders, const double scale, const MatrixI& idBlocks) const
   {
      assert(idBlocks.rows() > 0);
      assert(this->mPosN + this->mNegN <= in.rows());

      this->mTmp.topRows(this->mPosN).setZero();
      this->mTmp.bottomRows(this->mNegN).setZero();

      bool isComplex;
      double sgn;
      Array fp(this->mPosN);
      Array fn(this->mNegN);
      ArrayZ sZp;
      ArrayZ sZn;
      for(auto& o: orders)
      {
         int fO = o.first;
         int sO = o.second;
         if(fO > 0)
         {
            if(fO%2 == 0)
            {
               sgn = std::pow(-1.0,(fO/2)%2);
               isComplex = false;
            } else
            {
               sgn = std::pow(-1.0,((fO-1)/2)%2);
               isComplex = true;
            }
            fp = (sgn*(scale*this->positiveK()).array().pow(fO).matrix());
            fn = (sgn*(scale*this->negativeK()).array().pow(fO).matrix());
         } else
         {
            fp.setOnes(this->mPosN);
            fn.setOnes(this->mNegN);
            isComplex = false;
         }

         if(sO%2 == 0)
         {
            sgn = std::pow(-1.0,(sO/2)%2);
         } else
         {
            sgn = std::pow(-1.0,((sO-1)/2)%2);
         }
         isComplex = isComplex ^ (sO%2 == 1);

         if(isComplex)
         {
            sZp.resize(fp.size());
            sZp.real().setZero();
            sZn.resize(fn.size());
            sZn.real().setZero();
         }

         int start = 0;
         int negRow = this->mTmp.rows() - this->mNegN;
         for(int i = 0; i < idBlocks.rows(); ++i)
         {
            double k = sgn*std::pow(scale*idBlocks(i,0),sO);
            if((fO + sO)%2 == 1)
            {
               sZp.imag() = k*fp;
               sZn.imag() = k*fn;
            } else if(fO%2 == 1)
            {
               k = -k;
            }

            // Split positive and negative frequencies
            if(isComplex)
            {
               this->mTmp.block(0, start, this->mPosN, idBlocks(i,1)) += sZp.asDiagonal()*in.block(0, start, this->mPosN, idBlocks(i,1));
               this->mTmp.block(negRow, start, this->mNegN, idBlocks(i,1)) += sZn.asDiagonal()*in.block(this->mPosN, start, this->mNegN, idBlocks(i,1));
            } else
            {
               this->mTmp.block(0, start, this->mPosN, idBlocks(i,1)) += (k*fp.array()).matrix().asDiagonal()*in.block(0, start, this->mPosN, idBlocks(i,1));
               this->mTmp.block(negRow, start, this->mNegN, idBlocks(i,1)) += (k*fn.array()).matrix().asDiagonal()*in.block(this->mPosN, start, this->mNegN, idBlocks(i,1));
            }

            // Increment block counter
            start += idBlocks(i,1);
         }
      }

      this->forceConjugate(this->mTmp);
      this->applyPadding(this->mTmp);
   }

   void ComplexProjector::forceConjugate(MatrixZ& rData) const
   {
      int endN = rData.rows();

      // Loop over mean blocks
      for(auto it = this->mMeanBlocks.cbegin(); it != this->mMeanBlocks.cend(); ++it)
      {
         // Copy complex conjugate into negative frequency part
         for(int i = 1; i < this->mPosN; i++)
         {
            rData.block(endN - i, it->first, 1, it->second) = rData.block(i, it->first, 1, it->second).conjugate();
         }
      }
   }

   void ComplexProjector::applyPadding(MatrixZ& rData) const
   {
      rData.block(this->mPosN, 0, this->mPadSize, rData.cols()).setZero();
   }

}
}
}
}
}
