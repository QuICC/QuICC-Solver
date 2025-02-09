/**
 * @file ComplexIntegrator.cu
 * @brief Source of the interface for a generic cuFFT based complex integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/ComplexIntegrator.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   ComplexIntegrator::ComplexIntegrator()
      : IComplexBackend(1, 1), mNextDiff2DId(0), mOutMap(NULL,0,0)
   {
   }

   ComplexIntegrator::~ComplexIntegrator()
   {
   }

   void ComplexIntegrator::init(const SetupType& setup) const
   {
      //Initialize parent
      IComplexBackend::init(setup);

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
         cudaMalloc((void**)&(this->mcuFwd.back()), sizeof(cufftDoubleComplex)*this->mFwdSize*this->mBlockSize);

         // Create CUDA stream
         cudaStream_t stream;
         cudaStreamCreate(&stream);
         this->mStreams.push_back(stream);

         // Create cuFFT plan
         cufftHandle plan;
         // Create the complex to real plan
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_Z2Z, this->mBlockSize) != CUFFT_SUCCESS)
         {
            throw std::logic_error("CUFFT Error: Unable to create plan");
         }
         this->mPlans.push_back(plan);

         // Assign stream
         cufftSetStream(this->mPlans.back(), this->mStreams.back());
      }
   }

   void ComplexIntegrator::io(MatrixZ& rOut, const MatrixZ& in) const
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

   void ComplexIntegrator::io(MHDComplex* out, const MHDComplex* in) const
   {
      IComplexBackend::io(out, in);

      new (&this->mOutMap) Eigen::Map<MatrixZ>(this->mpOut, this->mBwdSize, this->mBlockSize);
   }

   void ComplexIntegrator::output(MatrixZ& rOut) const
   {
      rOut.topRows(this->mPosN) = this->mFftScaling*this->mOutMap.topRows(this->mPosN);
      rOut.bottomRows(this->mNegN) = this->mFftScaling*this->mOutMap.bottomRows(this->mNegN);
   }

   void ComplexIntegrator::outputDiff(MatrixZ& rOut, const int order, const double scale) const
   {
      // Odd order is complex
      if(order%2 == 1)
      {
         MHDComplex sgn = std::pow(-1.0,((order-1)/2)%2)*this->mFftScaling*Math::cI;
         rOut.topRows(this->mPosN) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*this->mOutMap.topRows(this->mPosN);
         rOut.bottomRows(this->mNegN) = (sgn*(scale*this->negativeK()).array().pow(order).matrix()).asDiagonal()*this->mOutMap.bottomRows(this->mNegN);
      } else
      {
         double sgn = std::pow(-1.0,(order/2)%2)*this->mFftScaling;
         rOut.topRows(this->mPosN) = (sgn*(scale*this->positiveK()).array().pow(order).matrix()).asDiagonal()*this->mOutMap.topRows(this->mPosN);
         rOut.bottomRows(this->mNegN) = (sgn*(scale*this->negativeK()).array().pow(order).matrix()).asDiagonal()*this->mOutMap.bottomRows(this->mNegN);
      }
   }

   int ComplexIntegrator::addDiff2DOp(const MatrixI& idBlocks) const
   {
      this->mDiff2DOp.insert(std::pair<int,Diff2DOpType>(this->mNextDiff2DId, Diff2DOpType()));
      int id = this->mDiff2DOp.size()-1;
      this->mNextDiff2DId = this->mDiff2DOp.size();
      std::get<0>(this->mDiff2DOp.at(id)).setZero(this->mPosN, idBlocks.rows());
      std::get<1>(this->mDiff2DOp.at(id)).setZero(this->mNegN, idBlocks.rows());
      std::get<2>(this->mDiff2DOp.at(id)) = idBlocks;

      return id;
   }

   void ComplexIntegrator::destroyDiff2D(const int id) const
   {
      this->mDiff2DOp.erase(id);
      this->mNextDiff2DId = id;
   }

   int ComplexIntegrator::computeDiff2D(const std::vector<std::pair<int,int> >& orders, const double scale, const MatrixI& idBlocks, const bool isInverse) const
   {
      assert(idBlocks.rows() > 0);

      bool isComplex;
      double sgn;
      Array fp(this->mPosN);
      Array fn(this->mNegN);
      int id = this->addDiff2DOp(idBlocks);
      MatrixZ& opP = std::get<0>(this->mDiff2DOp.at(id));
      MatrixZ& opN = std::get<1>(this->mDiff2DOp.at(id));
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

         int start = 0;
         for(int i = 0; i < idBlocks.rows(); ++i)
         {
            double k = sgn*std::pow(scale*idBlocks(i,0),sO);
            if((fO + sO)%2 == 1)
            {
               opP.col(i).imag() += k*fp;
               opN.col(i).imag() += k*fn;
            } else
            {
               if(fO%2 == 1)
               {
                  k = -k;
               }
               opP.col(i).real() += k*fp;
               opN.col(i).real() += k*fn;
            }

            // Increment block counter
            start += idBlocks(i,1);
         }
      }

      if(isInverse)
      {
         opP = opP.cwiseInverse();
         opN = opN.cwiseInverse();
      }

      return id;
   }

   int ComplexIntegrator::multDiff2D(const int idA, const int idB) const
   {
      int id = this->addDiff2DOp(std::get<2>(this->mDiff2DOp.at(idA)));
      std::get<0>(this->mDiff2DOp.at(id)) = std::get<0>(this->mDiff2DOp.at(idA)).array()*std::get<0>(this->mDiff2DOp.at(idB)).array();
      std::get<1>(this->mDiff2DOp.at(id)) = std::get<1>(this->mDiff2DOp.at(idA)).array()*std::get<1>(this->mDiff2DOp.at(idB)).array();

      return id;
   }

   void ComplexIntegrator::applyDiff2D(MatrixZ& rOut, const int id) const
   {
      int negRow = this->mTmp.rows() - this->mNegN;
      MatrixZ& opP = std::get<0>(this->mDiff2DOp.at(id));
      MatrixZ& opN = std::get<1>(this->mDiff2DOp.at(id));
      MatrixI& idBlocks = std::get<2>(this->mDiff2DOp.at(id));

      int start = 0;
      for(int i = 0; i < idBlocks.rows(); ++i)
      {
         // Split positive and negative frequencies
         this->mOutMap.block(0, start, this->mPosN, idBlocks(i,1)) = opP.col(i).asDiagonal()*this->mOutMap.block(0, start, this->mPosN, idBlocks(i,1));
         this->mOutMap.block(negRow, start, this->mNegN, idBlocks(i,1)) = opN.col(i).asDiagonal()*this->mOutMap.block(negRow, start, this->mNegN, idBlocks(i,1));

         // Increment block counter
         start += idBlocks(i,1);
      }

      this->output(rOut);
   }

   void ComplexIntegrator::outputZeroMean(MatrixZ& rOut) const
   {
      // Zero the mean
      for(auto it = this->mMeanBlocks.cbegin(); it != this->mMeanBlocks.cend(); ++it)
      {
         this->mOutMap.block(1, it->first, this->mOutMap.rows()-1, it->second).setZero();
      }

      this->output(rOut);
   }

   void ComplexIntegrator::outputMean(MatrixZ& rOut) const
   {
      rOut.setZero();

      // Set the mean
      for(auto it = this->mMeanBlocks.cbegin(); it != this->mMeanBlocks.cend(); ++it)
      {
         rOut.block(0, it->first, 1, it->second) = this->mFftScaling*this->mOutMap.block(0, it->first, 1, it->second);
      }
   }

   void ComplexIntegrator::applyFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNStreams;
         cudaMemcpyAsync(this->mcuFwd.at(sid), this->mpIn+fshift, sizeof(cufftDoubleComplex)*this->mFwdSize*this->mBlockSize,cudaMemcpyHostToDevice, this->mStreams.at(sid));
         fshift += this->mFwdSize*this->mBlockSize;
         cufftExecZ2Z(this->mPlans.at(sid), this->mcuFwd.at(sid), this->mcuBwd.at(sid), CUFFT_FORWARD);
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

   void ComplexIntegrator::extractMean() const
   {
      for(auto it = this->mMeanBlocks.cbegin(); it != this->mMeanBlocks.cend(); ++it)
      {
         this->mTmpMean.push_back(this->mOutMap.block(0, it->first, 1, it->second).transpose());
      }
   }

   void ComplexIntegrator::setMean(MatrixZ& rOut, const double scale) const
   {
      double f = scale*this->mFftScaling;
      for(size_t j = 0; j < this->mTmpMean.size(); ++j)
      {
         rOut.block(0, this->mMeanBlocks.at(j).first, 1, this->mMeanBlocks.at(j).second) = f*this->mTmpMean.at(j).transpose();
      }
      this->mTmpMean.clear();
   }

}
}
}
}
}
