/**
 * @file WorlandIntegrator.cu
 * @brief Source of the interface for a generic cuFFT based Worland integrator
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
#include "QuICC/Transform/Fft/Backend/CuFft/WorlandIntegrator.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuDctTools.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuMatrix.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/CheckCuda.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   WorlandIntegrator::WorlandIntegrator()
      : IWorlandBackend(1,1)
   {
   }

   WorlandIntegrator::~WorlandIntegrator()
   {
   }

   void WorlandIntegrator::init(const SetupType& setup, const int lshift, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      //Initialize parent
      IWorlandBackend::init(setup, lshift, lshiftOnlyParity, alwaysZeroNegative);
      this->mpcuWorkTmp = &this->mcuOutTmp;

      // Set internal spectral resolution
      this->setWSize(0);
      assert(this->mWSize <= this->mBwdSize);

      // Set transform scaling
      this->mFftScaling = std::sqrt(Math::PI)/static_cast<double>(this->mFwdSize);

      // Initialize main temporary storage
      int blockSize = std::max(this->mEBlockSize,this->mOBlockSize);
      // Input temporary storage
      this->initStorage(this->mFwdSize, blockSize, 1, this->mcuInTmp);
      // Output temporary storage
      this->initStorage(this->mBwdSize, blockSize, 2, this->mcuOutTmp);

      this->mPinnedIn.resize(2*std::max(this->mFwdSize,this->mBwdSize),(this->mEBlockSize+this->mOBlockSize));
      this->mPinnedOut.resize(2*std::max(this->mFwdSize,this->mBwdSize),(this->mEBlockSize+this->mOBlockSize));

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
      int cuBatch = this->mEBlockSize/this->mNBatches + (this->mEBlockSize%this->mNBatches > 0);
      for(int i = 0; i < this->mNBatches; i++)
      {
         // Initialise GPU storage
         this->mcuFwd.push_back(0);
         this->mcuBwd.push_back(0);
         this->mcuWork.push_back(0);
         CheckCuda(cudaMalloc((void**)&(this->mcuBwd.back()), sizeof(cufftDoubleComplex)*cuBwdSize*cuBlockSize), __LINE__);
         CheckCuda(cudaMalloc((void**)&(this->mcuFwd.back()), sizeof(cufftDoubleReal)*this->mFwdSize*cuBlockSize), __LINE__);
         CheckCuda(cudaMalloc((void**)&(this->mcuWork.back()), sizeof(cufftDoubleReal)*this->mFwdSize*cuBlockSize), __LINE__);

         // Create cuFFT plan
         cufftHandle plan;
         // Create the complex to real plan
         if(cufftPlanMany(&plan, 1, fftSize, NULL, 1, 0, NULL, 1, 0, CUFFT_D2Z, cuBatch) != CUFFT_SUCCESS)
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
      cuBatch = this->mOBlockSize/this->mNBatches + (this->mOBlockSize%this->mNBatches > 0);
      for(int i = 0; i < this->mNBatches; i++)
      {
         // Initialise GPU storage
         this->mcuWorkZ.push_back(0);
         CheckCuda(cudaMalloc((void**)&(this->mcuWorkZ.back()), sizeof(cufftDoubleComplex)*fwdSize*cuBlockSize), __LINE__);

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

      this->io(this->mcuOutTmp.at(0).data(), this->mcuInTmp.at(0).data());
   }

   void WorlandIntegrator::input(const Matrix& in, const bool isEven) const
   {
      GpuMatrix& cuIn = this->mcuInTmp.at(0);
      int start = 0;
      for(auto loc: *this->pLoc(isEven))
      {
         cublasSetMatrix(in.rows(), std::get<2>(loc), sizeof(double), in.data()+std::get<1>(loc)*in.rows(), in.rows(), cuIn.data()+start*cuIn.rows(), cuIn.rows());
         start += std::get<2>(loc);
      }
   }

#if 0
   void WorlandIntegrator::input(const MatrixZ& in, const bool isEven, const bool useReal) const
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

      for(auto loc: *this->pLoc(isEven))
      {
         for(int k = 0; k < std::get<2>(loc); k++)
         {
            cublasSetVector(in.rows(), sizeof(double), reinterpret_cast<const double*>(in.data())+2*(std::get<1>(loc)+k)*in.rows()+zShift, 2, cuIn.data()+(start+k)*cuIn.rows(), 1);
         }
         start += std::get<2>(loc);
      }
   }
#else
#if 0
   void WorlandIntegrator::input(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      GpuMatrix& cuIn = this->mcuInTmp.at(0);

      int zShift;
      if(useReal)
      {
         zShift = 0;
      } else
      {
         zShift = 1;
      }

      int start = 0;
      for(auto loc: *this->pLoc(isEven))
      {
         int cols = std::get<2>(loc);
         int iStart = 2*std::get<1>(loc)*in.rows();
         memcpy(this->mPinnedIn.data(), reinterpret_cast<const double*>(in.data())+iStart, 2*sizeof(double)*in.rows()*cols);
         CheckCuda(cudaMemcpy2D(cuIn.data() + start*cuIn.rows(), sizeof(double), this->mPinnedIn.data()+zShift, 2*sizeof(double), sizeof(double), in.rows()*cols, cudaMemcpyHostToDevice), __LINE__);
         start += std::get<2>(loc);
      }

   }
#else
   void WorlandIntegrator::input(const MatrixZ& in, const bool isEven, const bool useReal) const
   {
      GpuMatrix& cuIn = this->mcuInTmp.at(0);

      int zShift;
      if(useReal)
      {
         zShift = 0;
      } else
      {
         zShift = 1;
      }

      if(useReal)
      {
         int start = 0;
         for(auto loc: *this->pLoc(isEven))
         {
            int cols = std::get<2>(loc);
            int iStart = 2*std::get<1>(loc)*in.rows();
            memcpy(this->mPinnedIn.data()+2*start*in.rows(), reinterpret_cast<const double*>(in.data())+iStart, 2*sizeof(double)*in.rows()*cols);
            start += cols;
         }
      }

      CheckCuda(cudaMemcpy2D(cuIn.data(), sizeof(double), this->mPinnedIn.data()+zShift, 2*sizeof(double), sizeof(double), cuIn.rows()*this->blockSize(isEven), cudaMemcpyHostToDevice), __LINE__);

   }
#endif
#endif

   void WorlandIntegrator::applyEvenFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      dim3 kblocks(16,16);

      double *sn, *cs;
      CheckCuda(cudaMalloc((void**)&sn, sizeof(double)*this->mBwdSize), __LINE__);
      CheckCuda(cudaMalloc((void**)&cs, sizeof(double)*this->mBwdSize), __LINE__);
      Gpu::buildDCTPhases<<<1,256>>>(sn, cs, 2*this->mBwdSize, this->mBwdSize);

      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNBatches;
         int blocks = this->mcuEBlockSize.at(i);
         Gpu::buildDCTInput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuWork.at(sid),this->mpIn+fshift, this->mFwdSize, blocks);
         fshift += this->mFwdSize*blocks;
         cufftExecD2Z(this->mPlans.at(sid), this->mcuWork.at(sid), this->mcuBwd.at(sid));
         Gpu::buildDCTOutput<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mpOut+bshift,this->mcuBwd.at(sid), this->mBwdSize, blocks, sn, cs);
         bshift += this->mBwdSize*blocks;
      }
      CheckCuda(cudaFree(sn), __LINE__);
      CheckCuda(cudaFree(cs), __LINE__);

      for(int i = 0; i < this->mNBatches; i++)
      {
        this->synchronize(i);
      }
   }

   void WorlandIntegrator::applyOddFft() const
   {
      int fshift = 0;
      int bshift = 0;
      int sid = 0;
      dim3 kblocks(16,16);

      double *sn, *cs;
      CheckCuda(cudaMalloc((void**)&sn, sizeof(double)*this->mBwdSize), __LINE__);
      CheckCuda(cudaMalloc((void**)&cs, sizeof(double)*this->mBwdSize), __LINE__);
      Gpu::buildDCTPhases<<<1,256>>>(sn, cs, this->mBwdSize, this->mBwdSize);

      double *sn4, *cs4;
      CheckCuda(cudaMalloc((void**)&sn4, sizeof(double)*this->mBwdSize), __LINE__);
      CheckCuda(cudaMalloc((void**)&cs4, sizeof(double)*this->mBwdSize), __LINE__);
      Gpu::buildDCT4Phases<<<1,256>>>(sn4, cs4, 2*this->mBwdSize, this->mBwdSize);

      for(int i = 0; i < this->mNBatches; i++)
      {
         sid = i % this->mNBatches;
         int blocks = this->mcuOBlockSize.at(i);
         Gpu::buildDCT4Input<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mcuBwd.at(sid),this->mpIn+bshift, this->mBwdSize, blocks, sn, cs);
         bshift += this->mBwdSize*blocks;
         cufftExecZ2Z(this->mPlans.at(sid+this->mNBatches), this->mcuBwd.at(sid), this->mcuWorkZ.at(sid), CUFFT_FORWARD);
         Gpu::buildDCT4Output<<<1,kblocks,0,this->mStreams.at(sid)>>>(this->mpOut+fshift,this->mcuWorkZ.at(sid), this->mFwdSize, blocks, 1.0, sn4, cs4);
         fshift += this->mFwdSize*blocks;
      }
      CheckCuda(cudaFree(sn), __LINE__);
      CheckCuda(cudaFree(cs), __LINE__);
      CheckCuda(cudaFree(sn4), __LINE__);
      CheckCuda(cudaFree(cs4), __LINE__);

      for(int i = 0; i < this->mNBatches; i++)
      {
        this->synchronize(i);
      }
   }

   void WorlandIntegrator::applyFft() const
   {
     if(this->mPlanIsEven)
     {
       this->applyEvenFft();
     } else
     {
       this->applyOddFft();
     }
   }

   void WorlandIntegrator::partialForwardWorland(const int l, const int i0, const int start, const int cols, const GpuMatrix& J, double* out, const int oSize, const double normV, const double normM) const
   {
      int sId;
      for(int i = i0; i < l/2; ++i)
      {
         sId = i%this->mNStreams;
         this->mcuV.reshape(GpuMatrix::KEEP_SIZE, this->lSize(i,l));
         this->mcuM.reshape(GpuMatrix::KEEP_SIZE, this->lSize(i,l));
         this->buildShiftPair(sId, this->mcuV, this->mcuM, i, J, normV, normM);
         this->synchronize(sId);
         this->applyPair(sId, out, oSize, start, cols, this->mcuV, this->mcuM);
      }
   }

   void WorlandIntegrator::forwardWorland(const bool isEven, const int id) const
   {
      // Reset current l
      this->resetLocations(isEven, id);

      const GpuMatrix& J = this->J(isEven);
      GpuMatrix& cuOut = this->cuWorkTmp(id);

      const double normV = 1.0;
      const double normM = 1.0;
      int start = 0;
      int i0 = 0;
      int cols = cuOut.cols();
      double scale;
      const double zero = 0;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         cols = cuOut.cols() - start;
         if(l < 2)
         {
            scale = this->mFftScaling;
            double *pO = cuOut.data()+start*cuOut.rows();
            int ldo = cuOut.rows();
            CheckCuda(cublasDgeam(this->mHBlas, CUBLAS_OP_N, CUBLAS_OP_N, this->lSize(l), std::get<2>(loc), &scale, pO, ldo, &zero, pO, ldo, pO, ldo), __LINE__);
            if(l == 0)
            {
               scale = 1.0/std::sqrt(2.0);
               cublasDscal(this->mHBlas, std::get<2>(loc), &scale, cuOut.data()+start*cuOut.rows(), cuOut.rows());
            }
         } else
         {
            if(i0 == 0)
            {
               scale = this->mFftScaling;
               double *pO = cuOut.data()+start*cuOut.rows();
               int ldo = cuOut.rows();
               CheckCuda(cublasDgeam(this->mHBlas, CUBLAS_OP_N, CUBLAS_OP_N, this->lSize(0,l), cols, &scale, pO, ldo, &zero, pO, ldo, pO, ldo), __LINE__);
            }
            this->partialForwardWorland(l, i0, start, cols, J, cuOut.data(), cuOut.rows(), normV, normM);
            i0 = l/2;
         }
         std::get<3>(loc) = this->lSize(l);
         start += std::get<2>(loc);
      }
      this->synchronize();
   }

   void WorlandIntegrator::lowerBeta(const double alpha, const bool isEven, const int id, const double norm) const
   {
      GpuMatrix& cuOut = this->cuWorkTmp(id);

      int start = 0;
      double scale;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         this->mcuV.reshape(GpuMatrix::KEEP_SIZE, this->lSize(l));
         this->buildShiftV(this->mcuV, alpha, l-0.5, norm);
         this->applyTriSolve(cuOut.data(), cuOut.rows(), start, cols, this->mcuV);
         if(l == 1)
         {
            scale = 1.0/std::sqrt(2.0);
            cublasDscal(this->mHBlas, std::get<2>(loc), &scale, cuOut.data()+start*cuOut.rows(), cuOut.rows());
         }
         std::get<4>(loc)--;
         start += cols;
      }
   }

   void WorlandIntegrator::raiseBeta(const double alpha, const bool isEven, const int id, const double norm) const
   {
      GpuMatrix& cuOut = this->cuWorkTmp(id);

      int start = 0;
      double scale;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l < 0)
         {
            scale = 0.0;
            cublasDscal(this->mHBlas, cuOut.rows()*cols, &scale, cuOut.data()+start*cuOut.rows(), 1);
         } else
         {
            this->mcuV.reshape(GpuMatrix::KEEP_SIZE, this->lSize(l));
            this->buildShiftV(this->mcuV, alpha, l+0.5, norm);
            if(l==0)
            {
               scale = std::sqrt(2.0);
               cublasDscal(this->mHBlas, 1, &scale, this->mcuV.data()+1, 1);
            }
            this->applyTriProduct(cuOut.data(), cuOut.rows(), start, cols, this->mcuV);
            std::get<4>(loc)++;
         }
         std::get<3>(loc)--;
         start += cols;
      }
   }

   void WorlandIntegrator::lowerR2Beta(const double alpha, const bool isEven, const int id, const double norm) const
   {
      GpuMatrix& cuOut = this->cuWorkTmp(id);
      int start = 0;
      double scale;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l < 0)
         {
            scale = 0.0;
            cublasDscal(this->mHBlas, cuOut.rows()*cols, &scale, cuOut.data()+start*cuOut.rows(), 1);
         } else
         {
            this->mcuM.reshape(GpuMatrix::KEEP_SIZE, this->lSize(l));
            this->buildShiftM(this->mcuM, alpha, l-0.5, norm);
            this->applyTriProduct(cuOut.data(), cuOut.rows(), start, cols, this->mcuM);
            std::get<4>(loc)--;
         }
         start += cols;
      }
   }

   void WorlandIntegrator::raiseR2Beta(const double alpha, const bool isEven, const int id, const double norm) const
   {
      throw std::logic_error("Raise r^2 beta operator has not been tested!");

      GpuMatrix& cuOut = this->cuWorkTmp(id);
      int start = 0;
      double scale;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         this->mcuM.reshape(GpuMatrix::KEEP_SIZE, this->lSize(l));
         this->buildShiftM(this->mcuM, alpha, l+0.5, norm);
         this->applyTriSolve(cuOut.data(), cuOut.rows(), start, cols, this->mcuM);
         if(l == 1)
         {
            scale = 1.0/std::sqrt(2.0);
            cublasDscal(this->mHBlas, std::get<2>(loc), &scale, cuOut.data()+start*cuOut.rows(), cuOut.rows());
         }
         std::get<4>(loc)++;
         start += cols;
      }
   }

   void WorlandIntegrator::lowerAlpha(const double alpha, const bool isEven, const int id, const double norm) const
   {
      GpuMatrix& cuOut = this->cuWorkTmp(id);
      int start = 0;
      double scale;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         this->mcuU.reshape(GpuMatrix::KEEP_SIZE, this->lSize(l));
         this->buildShiftU(this->mcuU, alpha, l-0.5, norm);
         this->applyTriSolve(cuOut.data(), cuOut.rows(), start, cols, this->mcuU);
         if(l == 1)
         {
            scale = 1.0/std::sqrt(2.0);
            cublasDscal(this->mHBlas, std::get<2>(loc), &scale, cuOut.data()+start*cuOut.rows(), cuOut.rows());
         }
         start += cols;
      }
   }

   void WorlandIntegrator::raiseAlpha(const double alpha, const bool isEven, const int id, const double norm) const
   {
      throw std::logic_error("Raise alpha operator has not been tested!");

      GpuMatrix& cuOut = this->cuWorkTmp(id);
      int start = 0;
      double scale;
      for(auto& loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         if(l < 0)
         {
            scale = 0.0;
            cublasDscal(this->mHBlas, cuOut.rows()*cols, &scale, cuOut.data()+start*cuOut.rows(), 1);
         } else
         {
            this->mcuU.reshape(GpuMatrix::KEEP_SIZE, this->lSize(l));
            this->buildShiftU(this->mcuU, alpha, l-0.5, norm);
            this->applyTriProduct(cuOut.data(), cuOut.rows(), start, cols, this->mcuU);
         }
         std::get<3>(loc)--;
         start += cols;
      }
   }

   void WorlandIntegrator::applyI2(const bool isEven, const int id) const
   {
      GpuMatrix& cuOut = this->cuWorkTmp(id);

      GpuMatrix cubd;
      if(this->pLoc(isEven,id))
      {
         int l = std::get<4>(*this->pLoc(isEven,id)->begin());
         int r = this->lSize(l);
         cubd = std::move(GpuMatrix(5,r));
      }

      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         int r = this->lSize(l);
         ::QuICC::SparseSM::Worland::I2 spasm(r, r, -0.5, -0.5, l);

         Matrix bd(5,r);
         spasm.buildOp(bd);
         cubd.reshape(GpuMatrix::KEEP_SIZE, r);
         cudaMemcpy(cubd.data(), bd.data(), sizeof(double)*bd.size(), cudaMemcpyHostToDevice);
         this->applyBandProduct(cuOut.data(), cuOut.rows(), start, cols, cubd);

         start += cols;
      }
   }

   void WorlandIntegrator::applyI4(const bool isEven, const int id) const
   {
      GpuMatrix& cuOut = this->cuWorkTmp(id);

      GpuMatrix cubd;
      if(this->pLoc(isEven,id))
      {
         int l = std::get<4>(*this->pLoc(isEven,id)->begin());
         int r = this->lSize(l);
         cubd = std::move(GpuMatrix(9,r));
      }

      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int l = std::get<4>(loc);
         int cols = std::get<2>(loc);
         int r = this->lSize(l);
         ::QuICC::SparseSM::Worland::I4 spasm(r, r, -0.5, -0.5, l);

         Matrix bd(9,r);
         spasm.buildOp(bd);
         cubd.reshape(GpuMatrix::KEEP_SIZE, r);
         cudaMemcpy(cubd.data(), bd.data(), sizeof(double)*bd.size(), cudaMemcpyHostToDevice);
         this->applyBandProduct(cuOut.data(), cuOut.rows(), start, cols, cubd);

         start += cols;
      }
   }

   void WorlandIntegrator::output(Matrix& rOut, const bool isEven) const
   {
      GpuMatrix& cuOut = this->mcuOutTmp.at(0);
      int start = 0;
      for(auto loc: *this->pLoc(isEven))
      {
         int cols = std::get<2>(loc);
         for(int k = 0; k < cols; k++)
         {
            cublasGetVector(this->mSpecSize, sizeof(double), cuOut.data()+(start+k)*cuOut.rows(), 1, rOut.data()+(std::get<1>(loc)+k)*rOut.rows(), 1);
         }
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

#if 0
   void WorlandIntegrator::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      GpuMatrix& cuOut = this->mcuOutTmp.at(0);
      int start = 0;
      if(useReal)
      {
         for(auto loc: *this->pLoc(isEven))
         {
            int cols = std::get<2>(loc);
            int zShift = 0;
            for(int k = 0; k < cols; k++)
            {
               cublasGetVector(this->mSpecSize, sizeof(double), cuOut.data()+(start+k)*cuOut.rows(), 1, reinterpret_cast<double*>(rOut.data())+2*(std::get<1>(loc)+k)*rOut.rows()+zShift, 2);
            }
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
            int zShift = 1;
            for(int k = 0; k < cols; k++)
            {
               cublasGetVector(this->mSpecSize, sizeof(double), cuOut.data()+(start+k)*cuOut.rows(), 1, reinterpret_cast<double*>(rOut.data())+2*(std::get<1>(loc)+k)*rOut.rows()+zShift, 2);
            }
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
#else
#if 0
   void WorlandIntegrator::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      GpuMatrix& cuOut = this->mcuOutTmp.at(0);
      int start = 0;
      const double one = 1.0;
      const double zero = 0.0;

      assert(rOut.rows() == this->mSpecSize);

      GpuMatrix& cuWork = this->cuWorkBuffer();
      const double *pO = cuOut.data();
      int ldo = cuOut.rows();
      CheckCuda(cublasDgeam(this->mHBlas, CUBLAS_OP_N, CUBLAS_OP_N, this->mSpecSize, this->blockSize(isEven), &one, pO, ldo, &zero, pO, ldo, cuWork.data(), this->mSpecSize), __LINE__);

      int zShift;
      if(useReal)
      {
         zShift = 0;
      } else
      {
         zShift = 1;
      }

      for(auto loc: *this->pLoc(isEven))
      {
         int cols = std::get<2>(loc);
         int oStart = 2*std::get<1>(loc)*rOut.rows();
         CheckCuda(cudaMemcpy2D(this->mPinnedOut.data()+oStart+zShift, 2*sizeof(double), cuWork.data() + start*this->mSpecSize, sizeof(double), sizeof(double), this->mSpecSize*cols, cudaMemcpyDeviceToHost), __LINE__);
         start += cols;
      }

      if(!useReal && !isEven)
      {
         memcpy(reinterpret_cast<double*>(rOut.data()), this->mPinnedOut.data(), 2*sizeof(double)*this->mSpecSize*rOut.cols());

         // Zero unused values
         for(auto loc: this->mZLoc)
         {
            int s = std::get<1>(loc);
            int cols = std::get<2>(loc);
            rOut.block(0, s, this->mSpecSize, cols).setZero();
         }
      }

   }
#else
   void WorlandIntegrator::output(MatrixZ& rOut, const bool isEven, const bool useReal) const
   {
      GpuMatrix& cuOut = this->mcuOutTmp.at(0);
      int start = 0;
      const double one = 1.0;
      const double zero = 0.0;

      assert(rOut.rows() == this->mSpecSize);

      GpuMatrix& cuWork = this->cuWorkBuffer();
      const double *pO = cuOut.data();
      int ldo = cuOut.rows();
      CheckCuda(cublasDgeam(this->mHBlas, CUBLAS_OP_N, CUBLAS_OP_N, this->mSpecSize, this->blockSize(isEven), &one, pO, ldo, &zero, pO, ldo, cuWork.data(), this->mSpecSize), __LINE__);

      int zShift;
      if(useReal)
      {
         zShift = 0;
      } else
      {
         zShift = 1;
      }
      CheckCuda(cudaMemcpy2D(this->mPinnedOut.data()+zShift, 2*sizeof(double), cuWork.data(), sizeof(double), sizeof(double), this->mSpecSize*this->blockSize(isEven), cudaMemcpyDeviceToHost), __LINE__);

      if(!useReal)
      {
         for(auto loc: *this->pLoc(isEven))
         {
            int cols = std::get<2>(loc);
            int oStart = 2*std::get<1>(loc)*rOut.rows();
            memcpy(reinterpret_cast<double*>(rOut.data())+oStart, this->mPinnedOut.data()+2*start*rOut.rows(), 2*sizeof(double)*rOut.rows()*cols);
            start += cols;
         }
      }

      if(!useReal && !isEven)
      {
         // Zero unused values
         for(auto loc: this->mZLoc)
         {
            int s = std::get<1>(loc);
            int cols = std::get<2>(loc);
            rOut.block(0, s, this->mSpecSize, cols).setZero();
         }
      }

   }
#endif
#endif

}
}
}
}
}
