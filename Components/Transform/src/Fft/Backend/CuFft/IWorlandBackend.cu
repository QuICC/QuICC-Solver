/**
 * @file IWorlandBackend.cu
 * @brief Source of the interface for a generic cuFFT based Worland integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/CuFft/CheckCuda.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/IWorlandBackend.hpp"

// Project includes
//

// Banded triangular product implementation
//#define QUICC_TRIPRODUCT_TBMV
#define QUICC_TRIPRODUCT_SPMM

// Banded triangular solve implementation
//#define QUICC_TRISOLVE_TBSV
//#define QUICC_TRISOLVE_CSRSV2
#define QUICC_TRISOLVE_CSRSM2

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   const int IWorlandBackend::MIN_BATCH_BLOCKSIZE = 10;

   const double IWorlandBackend::UPPER_BANDED = 424242.424242;

   const double IWorlandBackend::LOWER_BANDED = 434343.434343;

   const double IWorlandBackend::LOWER_BANDED_PADDED = -IWorlandBackend::LOWER_BANDED;

   IWorlandBackend::IWorlandBackend()
      : IWorlandBackend(1,1)
   {
   }

   IWorlandBackend::IWorlandBackend(const int nStreams, const int nBatches)
      : ICuFftBackend(nStreams, nBatches), mPlanIsEven(true), mFlipped(false), mSpecSize(-1), mFwdSize(-1), mBwdSize(-1), mWSize(-1), mEBlockSize(-1), mOBlockSize(-1)
   {
      this->mELoc.push_back(LocationVector());
      this->mOLoc.push_back(LocationVector());
   }

   IWorlandBackend::~IWorlandBackend()
   {
      for(auto& d: this->mcuFwd)
      {
         CheckCuda(cudaFree(d), __LINE__, false);
      }
      for(auto& d: this->mcuBwd)
      {
         CheckCuda(cudaFree(d), __LINE__, false);
      }
      for(auto& d: this->mcuWork)
      {
         CheckCuda(cudaFree(d), __LINE__, false);
      }
      for(auto& d: this->mcuWorkZ)
      {
         CheckCuda(cudaFree(d), __LINE__, false);
      }
      CheckCuda(cublasDestroy(this->mHBlas), __LINE__, false);
      CheckCuda(cusparseDestroy(this->mHSparse), __LINE__, false);
   }

   void IWorlandBackend::initStorage(const int rows, const int cols, const int n, std::vector<GpuMatrix>& s) const
   {
      assert(rows > 0);
      assert(cols > 0);
      assert(n > 0);

      // Input temporary storage
      s.reserve(n);
      s.push_back(GpuMatrix(rows, cols));
      CheckCuda(cudaMemset(s.back().data(), 0, sizeof(double)*s.back().size()), __LINE__);

      // Initialize additional temporary storage
      for(int i = 0; i < n-1; i++)
      {
         s.push_back(s.at(0));
      }
   }

   void IWorlandBackend::addStorage(const int inExtras, const int outExtras) const
   {
      assert(inExtras >= 0);
      assert(outExtras >= 0);

      // Initialize additional input temporary storage
      this->mcuInTmp.reserve(this->mcuInTmp.size() + inExtras);
      for(int i = 0; i < inExtras; i++)
      {
         this->mcuInTmp.push_back(this->mcuInTmp.at(0));
      }

      // Initialize additional put temporary storage
      this->mcuOutTmp.reserve(this->mcuOutTmp.size() + outExtras);
      for(int i = 0; i < outExtras; i++)
      {
         this->mcuOutTmp.push_back(this->mcuOutTmp.at(0));
      }

      // Initialize additional block locations
      this->mELoc.reserve(this->mELoc.size() + std::max(inExtras,outExtras));
      this->mOLoc.reserve(this->mOLoc.size() + std::max(inExtras,outExtras));
      for(int i = 0; i < std::max(inExtras,outExtras); i++)
      {
         this->mELoc.push_back(LocationVector());
         this->mOLoc.push_back(LocationVector());
      }
   }

   void IWorlandBackend::init(const SetupType& setup, const int lshift, const int extraN, const bool lshiftOnlyParity, const bool alwaysZeroNegative) const
   {
      this->mFlipped = (std::abs(lshift)%2 == 1);
      int lshift_parity = lshift;
      if(lshiftOnlyParity)
      {
         lshift_parity = 0;
      }
      int lshift_zero = 0;
      if(lshiftOnlyParity && alwaysZeroNegative)
      {
         lshift_zero = lshift;
      }

      this->mSpecSize = setup.specSize();
      this->mFwdSize = setup.fwdSize();
      this->mBwdSize = setup.bwdSize();

      if(this->mFwdSize%2 == 1)
      {
         throw std::logic_error("cuFFT backend only supports even resolutions");
      }

      // Compute even and odd block sizes
      this->mEBlockSize = 0;
      this->mOBlockSize = 0;
      int col = 0;
      for(int k = 0; k < setup.slowSize(); ++k)
      {
         int l = setup.slow(k);
         int loc_l = l+lshift_parity;
         LocationType loc = std::make_tuple(loc_l,col,setup.mult(k), -1, loc_l);
         if(this->inZFilter(l) || std::get<0>(loc) + lshift_zero < 0)
         {
            this->mZLoc.push_back(loc);
         } else
         {
            if((l + lshift) % 2 == 0)
            {
               this->mEBlockSize += std::get<2>(loc);
               this->mELoc.at(0).push_back(loc);
            } else
            {
               this->mOBlockSize += std::get<2>(loc);
               this->mOLoc.at(0).push_back(loc);
            }
         }
         col += std::get<2>(loc);
      }

      // Recompute maximum number of batches
      int blockSize = std::min(this->mEBlockSize, this->mOBlockSize);
      this->mNBatches = std::min(blockSize/IWorlandBackend::MIN_BATCH_BLOCKSIZE + 1, this->mNBatches);
      this->mNBatches = std::min(this->mNStreams, this->mNBatches);

      this->mcuEBlockSize.reserve(this->mNBatches);
      this->mcuOBlockSize.reserve(this->mNBatches);
      int rEB = this->mEBlockSize%this->mNBatches;
      int rOB = this->mOBlockSize%this->mNBatches;
      for(int i = 0; i < this->mNBatches; i++)
      {
         this->mcuEBlockSize.push_back(this->mEBlockSize/this->mNBatches + (rEB > 0));
         rEB--;
         this->mcuOBlockSize.push_back(this->mOBlockSize/this->mNBatches + (rOB > 0));
         rOB--;
      }

      // Create CUBLAS context
      CheckCuda(cublasCreate(&this->mHBlas), __LINE__);

      // Create CUSPARSE context
      CheckCuda(cusparseCreate(&this->mHSparse), __LINE__);

      this->mWExtra = extraN;
   }

   void IWorlandBackend::setZFilter(const std::set<int>& filter) const
   {
      this->mZFilter = filter;
   }

   bool IWorlandBackend::isPhysEven(const bool isSpecEven) const
   {
      return (isSpecEven ^ this->mFlipped);
   }

   bool IWorlandBackend::inZFilter(const int l) const
   {
      return (this->mZFilter.count(l) == 1);
   }

   void IWorlandBackend::setWSize() const
   {
      // Triangular truncation requirements: WSize = n_lmax_modes + lmax/2
      this->mWSize = this->mSpecSize + this->mWExtra + (std::max(std::get<0>(*this->pLoc(true)->rbegin()),std::get<0>(*this->pLoc(false)->rbegin())) + shiftMaxL)/2;
   }

   IWorlandBackend::LocationVector* IWorlandBackend::pLoc(const bool isEven, const unsigned int id) const
   {
      LocationVector* ptr;
      if(this->isPhysEven(isEven))
      {
         assert(this->mELoc.size() > id);

         ptr = &this->mELoc.at(id);
      } else
      {
         assert(this->mOLoc.size() > id);

         ptr = &this->mOLoc.at(id);
      }

      return ptr;
   }

   void IWorlandBackend::resetLocations(const bool isEven, const int id) const
   {
      for(auto& loc: *this->pLoc(isEven,id))
      {
         std::get<4>(loc) = std::get<0>(loc);
      }
   }

   GpuMatrix& IWorlandBackend::cuWorkTmp(const unsigned int id) const
   {
      assert(this->mpcuWorkTmp);
      assert(this->mpcuWorkTmp->size() > id);

      return this->mpcuWorkTmp->at(id);
   }

   GpuMatrix& IWorlandBackend::cuWorkBuffer() const
   {
      assert(this->mpcuWorkTmp);
      assert(this->mpcuWorkTmp->size() > 1);

      return this->mpcuWorkTmp->back();
   }

   void IWorlandBackend::initJ() const
   {
      double alpha = -0.5;

      // Even modes
      if(this->mELoc.at(0).size() > 0)
      {
         // Max size is L + N + 1
         int lmax = std::get<0>(*this->mELoc.at(0).rbegin());
         this->mJEvenSize = this->lSize(lmax) + lmax  + 1;
         double beta = 0.5;
         this->mcuJEven = std::move(GpuMatrix(this->mJEvenSize, 4));
         Gpu::jacobiShiftMatrix<<<1,128,0,this->stream(0)>>>(this->mcuJEven.data(), this->mcuJEven.rows(), alpha, beta);
      }

      // Odd Modes
      if(this->mOLoc.at(0).size() > 0)
      {
         // Max size is L + N + 1
         int lmax = std::get<0>(*this->mOLoc.at(0).rbegin());
         this->mJOddSize = this->lSize(lmax) + lmax  + 1;
         double beta = 1.5;
         this->mcuJOdd = std::move(GpuMatrix(this->mJOddSize, 4));
         Gpu::jacobiShiftMatrix<<<1,128,0,this->stream(1%this->mNStreams)>>>(this->mcuJOdd.data(), this->mcuJOdd.rows(), alpha, beta);
      }
      this->synchronize();

      // Create GPU storage for operators
      this->mcuU.resize(2, this->mWSize+1);
      this->mcuV.resize(2, this->mWSize+1);
      this->mcuM.resize(2, this->mWSize+1);

      this->mcsrP = std::move(GpuCsrMatrix(this->mWSize+1, this->mWSize+1, 2*this->mWSize+1));
      this->mcsrS = std::move(GpuCsrMatrix(this->mWSize+1, this->mWSize+1, 2*this->mWSize+1));
      this->mcuBuffer = std::move(GpuMatrix(1024*1024, 1));
   }

   GpuMatrix& IWorlandBackend::J(const bool isEven) const
   {
      if(this->isPhysEven(isEven))
      {
         return this->mcuJEven;
      } else
      {
         return this->mcuJOdd;
      }
   }

   void IWorlandBackend::buildShiftU(GpuMatrix& U, const double alpha, const double beta, const double norm) const
   {
      GpuMatrix J(U.cols(), 4);
      Gpu::jacobiShiftMatrix<<<1,128>>>(J.data(), J.rows(), beta, alpha);
      Gpu::buildShiftU<<<1,128>>>(U.data(), 0, U.cols(), J.data(), J.rows(), norm);
   }

   void IWorlandBackend::buildShiftV(GpuMatrix& V, const double alpha, const double beta, const double norm) const
   {
      GpuMatrix J(V.cols(), 4);
      Gpu::jacobiShiftMatrix<<<1,128>>>(J.data(), J.rows(), alpha, beta);
      Gpu::buildShiftV<<<1,128>>>(V.data(), 0, V.cols(), J.data(), J.rows(), norm);
   }

   void IWorlandBackend::buildShiftM(GpuMatrix& M, const double alpha, const double beta, const double norm, const bool isSquare) const
   {
      GpuMatrix J(M.cols()+1, 4);
      Gpu::jacobiShiftMatrix<<<1,128>>>(J.data(), J.rows(), alpha, beta-1);
      Gpu::buildShiftM<<<1,128>>>(M.data(), 0, M.cols(), J.data(), J.rows(), norm, isSquare);
   }

   void IWorlandBackend::buildShiftPair(const int streamId, GpuMatrix& V, GpuMatrix& M, const int i, const GpuMatrix& J, const double normV, const double normM, const bool isSquare) const
   {
      Gpu::buildShiftV<<<1,128,0,this->stream(streamId)>>>(V.data(), i, V.cols(), J.data(), J.rows(), normV);

      Gpu::buildShiftM<<<1,128,0,this->stream(streamId)>>>(M.data(), i, M.cols(), J.data(), J.rows(), normM, isSquare);
   }

   void IWorlandBackend::setPlan(const bool isEven) const
   {
      this->mPlanIsEven = this->isPhysEven(isEven);
   }

   int IWorlandBackend::blockSize(const bool isEven) const
   {
      if(this->isPhysEven(isEven))
      {
         return this->mEBlockSize;
      } else
      {
         return this->mOBlockSize;
      }
   }

   void IWorlandBackend::applyTriProduct(double* out, const int oSize, const int start, const int cols, GpuMatrix& P) const
   {
      this->computeTriProduct(0, out, oSize, start, cols, P);
      this->synchronize(0);
   }

   void IWorlandBackend::applyBandProduct(double* out, const int oSize, const int start, const int cols, GpuMatrix& A) const
   {
      int r = A.cols();
      GpuMatrix& work = this->cuWorkBuffer();
      work.reshape(oSize, cols);
      CheckCuda(cublasDcopy(this->mHBlas, oSize*cols, out+start*oSize, 1, work.data(), 1), __LINE__);

      int kl, ku;
      this->getMatrixBands(A, kl, ku);

      cublasOperation_t trans = CUBLAS_OP_N;
      double alpha = 1.0;
      int lda = kl + ku + 1;
      int incx = 1;
      double beta = 0.0;
      for(int j = 0; j < cols; j++)
      {
         double *X = work.data() + j*oSize;
         double *Y = out + (start+j)*oSize;
         CheckCuda(cublasDgbmv(this->mHBlas, trans, r, r, kl, ku, &alpha, A.data(), lda, X, incx, &beta, Y, incx), __LINE__);
      }
      this->synchronize();
   }

   void IWorlandBackend::useStream(const int id) const
   {
      this->useStream(this->mHBlas, id);
      this->useStream(this->mHSparse, id);
   }

   void IWorlandBackend::useStream(cublasHandle_t h, const int id) const
   {
      if(id < 0)
      {
         CheckCuda(cublasSetStream(h, NULL), __LINE__);
      } else
      {
         CheckCuda(cublasSetStream(h, this->stream(id)), __LINE__);
      }
   }

   void IWorlandBackend::useStream(cusparseHandle_t h, const int id) const
   {
      if(id < 0)
      {
         CheckCuda(cusparseSetStream(h, NULL), __LINE__);
      } else
      {
         CheckCuda(cusparseSetStream(h, this->stream(id)), __LINE__);
      }
   }

   void IWorlandBackend::applyTriSolve(double* out, const int oSize, const int start, const int cols, GpuMatrix& S) const
   {
      this->computeTriSolve(0, out, oSize, start, cols, S);
      this->synchronize(0);
   }

   void IWorlandBackend::applyPair(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& P, GpuMatrix& S) const
   {
      this->computeTriProduct(streamId, out, oSize, start, cols, P);
      this->computeTriSolve(streamId, out, oSize, start, cols, S);
   }

   void IWorlandBackend::io(double* out, const double* in) const
   {
      this->mpOut = out;
      this->mpIn = in;
   }

   void IWorlandBackend::scaleC(const double c, const bool isEven, const unsigned int id) const
   {
      GpuMatrix& cuTmp = this->cuWorkTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         for(int k = 0; k < cols; k++)
         {
            CheckCuda(cublasDscal(this->mHBlas, rows, &c, cuTmp.data()+(start+k)*cuTmp.rows(), 1), __LINE__);
         }
         start += cols;
      }
      this->synchronize();
   }

   void IWorlandBackend::scaleALPY(const double a, const double y, const bool isEven, const int lshift, const unsigned int id) const
   {
      GpuMatrix& cuTmp = this->cuWorkTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven, id))
      {
         int l = std::get<4>(loc) + lshift;
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         double c = (a*l + y);
         for(int k = 0; k < cols; k++)
         {
            CheckCuda(cublasDscal(this->mHBlas, rows, &c, cuTmp.data()+(start+k)*cuTmp.rows(), 1), __LINE__);
         }
         start += cols;
      }
      this->synchronize();
   }

   void IWorlandBackend::scaleD(const bool isEven, const int lshift, const unsigned int id) const
   {
      const double alpha = -0.5;
      GpuMatrix& cuTmp = this->cuWorkTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven, id))
      {
         int l = std::get<4>(loc) + lshift;
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         for(int n = 0; n < rows; n++)
         {
            double c = 2.0*(l + n + 1.0)*std::sqrt((n+1.0)/(n + alpha + l + 1.5));
            CheckCuda(cublasDscal(this->mHBlas, cols, &c, cuTmp.data() + (start)*cuTmp.rows() + n, cuTmp.rows()), __LINE__);
         }
         start += cols;
      }
      this->synchronize();
   }

   void IWorlandBackend::scaleSphLaplA(const bool isEven, const int lshift, const unsigned int id) const
   {
      const double alpha = -0.5;
      GpuMatrix& cuTmp = this->cuWorkTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven, id))
      {
         int l = std::get<4>(loc) + lshift;
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         for(int n = 0; n < rows; n++)
         {
            double c = 4.0*(l + n + 2.0)*(l + n + 3.0)*std::sqrt((n+1.0)*(n+2.0)/((n + alpha + l + 2.5)*(n + alpha + l + 3.5)));
            CheckCuda(cublasDscal(this->mHBlas, cols, &c, cuTmp.data() + (start)*cuTmp.rows() + n, cuTmp.rows()), __LINE__);
         }
         start += cols;
      }
      this->synchronize();
   }

   void IWorlandBackend::scaleSphLaplB(const bool isEven, const int lshift, const unsigned int id) const
   {
      const double alpha = -0.5;
      GpuMatrix& cuTmp = this->cuWorkTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven, id))
      {
         int l = std::get<4>(loc) + lshift;
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         for(int n = 0; n < rows; n++)
         {
            double c = 2.0*(l + n + 1.0)*(2.0*l + 3.0)*std::sqrt((n+1.0)/(n + alpha + l + 1.5));
            CheckCuda(cublasDscal(this->mHBlas, cols, &c, cuTmp.data() + (start)*cuTmp.rows() + n, cuTmp.rows()), __LINE__);
         }
         start += cols;
      }
   }

   void IWorlandBackend::lshift(const unsigned int id, const int lshift, const bool isEven) const
   {
      for(auto& loc: *this->pLoc(isEven,id))
      {
         std::get<4>(loc) += lshift;
      }
   }

   void IWorlandBackend::nshift(const unsigned int id, const int nshift, const bool isEven) const
   {
      GpuMatrix& cuTmp = this->cuWorkTmp(id);
      int s = std::abs(nshift);
      GpuMatrix cuWork(cuTmp.rows(), 1);

      double zero = 0.0;
      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc) - s;
         if(nshift > 0)
         {
            for(int k = 0; k < cols; k++)
            {
               CheckCuda(cublasDcopy(this->mHBlas, rows, cuTmp.data() + (start+k)*cuTmp.rows(), 1, cuWork.data(), 1), __LINE__);
               CheckCuda(cublasDcopy(this->mHBlas, rows, cuWork.data(), 1, cuTmp.data() + (start+k)*cuTmp.rows() + s, 1), __LINE__);
               CheckCuda(cublasDscal(this->mHBlas, s, &zero, cuTmp.data() + (start+k)*cuTmp.rows(), 1), __LINE__);
            }
         } else
         {
            for(int k = 0; k < cols; k++)
            {
               CheckCuda(cublasDcopy(this->mHBlas, rows, cuTmp.data() + (start+k)*cuTmp.rows() + s, 1, cuWork.data(), 1), __LINE__);
               CheckCuda(cublasDcopy(this->mHBlas, rows, cuWork.data(), 1, cuTmp.data() + (start+k)*cuTmp.rows(), 1), __LINE__);
               CheckCuda(cublasDscal(this->mHBlas, s, &zero, cuTmp.data() + (start+k)*cuTmp.rows() + rows, 1), __LINE__);
            }
         }
         start += cols;
      }
   }

   void IWorlandBackend::copy(const int to, const int from, const int nshift, const bool isEven) const
   {
      assert(to != from);
      GpuMatrix& cuTmp = this->cuWorkTmp(from);
      GpuMatrix& cuExtra = this->cuWorkTmp(to);
      int s = std::abs(nshift);

      if(to > 0)
      {
         *this->pLoc(isEven,to) = *this->pLoc(isEven,from);
      }

      double zero = 0.0;
      int start = 0;
      for(auto loc: *this->pLoc(isEven,to))
      {
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc) - s;
         if(nshift > 0)
         {
            for(int k = 0; k < cols; k++)
            {
               CheckCuda(cublasDcopy(this->mHBlas, rows, cuTmp.data() + (start+k)*cuTmp.rows(), 1, cuExtra.data() + (start+k)*cuExtra.rows() + s, 1), __LINE__);
               CheckCuda(cublasDscal(this->mHBlas, s, &zero, cuExtra.data() + (start+k)*cuExtra.rows(), 1), __LINE__);
            }
         } else
         {
            for(int k = 0; k < cols; k++)
            {
               CheckCuda(cublasDcopy(this->mHBlas, rows, cuTmp.data() + (start+k)*cuTmp.rows() + s, 1, cuExtra.data() + (start+k)*cuExtra.rows(), 1), __LINE__);
               CheckCuda(cublasDscal(this->mHBlas, s, &zero, cuExtra.data() + (start+k)*cuExtra.rows()+ rows, 1), __LINE__);
            }
         }
         start += cols;
      }
      this->synchronize();
   }

   void IWorlandBackend::add(const int to, const int from, const int nshift, const bool isEven) const
   {
      GpuMatrix& cuTmp = this->cuWorkTmp(to);
      GpuMatrix& cuExtra = this->cuWorkTmp(from);
      int s = std::abs(nshift);

      double one = 1.0;
      int start = 0;
      for(auto loc: *this->pLoc(isEven,to))
      {
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc) - s;
         if(nshift > 0)
         {
            for(int k = 0; k < cols; k++)
            {
               CheckCuda(cublasDaxpy(this->mHBlas, rows, &one, cuExtra.data() + (start+k)*cuExtra.rows(), 1, cuTmp.data() + (start+k)*cuTmp.rows() + s, 1), __LINE__);
            }
         } else
         {
            for(int k = 0; k < cols; k++)
            {
               CheckCuda(cublasDaxpy(this->mHBlas, rows, &one, cuExtra.data() + (start+k)*cuExtra.rows() + s, 1, cuTmp.data() + (start+k)*cuTmp.rows(), 1), __LINE__);
            }
         }
         start += cols;
      }
      this->synchronize();
   }

   bool IWorlandBackend::getMatrixFillMode(const int streamId, const GpuMatrix& A, cublasFillMode_t& uplo) const
   {
      bool padded = false;

      double info[2];
      CheckCuda(cublasGetVectorAsync(2, sizeof(double), A.data(), A.size()-1, &info, 1, this->stream(streamId)), __LINE__);
      if(info[0] == UPPER_BANDED)
      {
         uplo = CUBLAS_FILL_MODE_UPPER;
      } else
      {
         assert(std::abs(info[1]) == LOWER_BANDED);
         uplo = CUBLAS_FILL_MODE_LOWER;
         if(info[1] == LOWER_BANDED_PADDED)
         {
            padded = true;
         }
      }

      return padded;
   }

   void IWorlandBackend::getMatrixBands(const GpuMatrix& A, int& kl, int& ku) const
   {
      double info[2];
      CheckCuda(cublasGetVector(2, sizeof(double), A.data(), A.size()-1, &info, 1), __LINE__);
      kl = static_cast<int>(info[1]);
      ku = static_cast<int>(info[0]);
   }

#ifdef QUICC_TRIPRODUCT_TBMV
   void IWorlandBackend::computeTriProduct(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& P) const
   {
      // Set stream
      this->useStream(streamId);

      int r = P.cols();
      const double zero = 0.0;

      cublasFillMode_t uplo;
      bool padded = this->getMatrixFillMode(P, uplo);
      if(padded)
      {
         CheckCuda(cublasDscal(this->mHBlas,cols,&zero,out+r+start*oSize,oSize), __LINE__);
      }

      cublasOperation_t trans = CUBLAS_OP_N;
      cublasDiagType_t diag = CUBLAS_DIAG_NON_UNIT;
      int k = 1;
      int lda = 2;
      int incx = 1;

      for(int j = 0; j < cols; j++)
      {
         double *X = out + (start+j)*oSize;
         CheckCuda(cublasDtbmv(this->mHBlas, uplo, trans, diag, r, k, P.data(), lda, X, incx), __LINE__);
      }

      // Set default stream
      this->useStream(-1);
   }
#endif //QUICC_TRIPRODUCT_TBMV

#ifdef QUICC_TRIPRODUCT_SPMM
   void IWorlandBackend::computeTriProduct(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& P) const
   {
      // Set stream
      this->useStream(streamId);

      int r = P.cols();
      const double zero = 0.0;
      const double one = 1.0;

      cublasFillMode_t P_uplo;
      bool padded = this->getMatrixFillMode(streamId, P, P_uplo);
      if(padded)
      {
         CheckCuda(cublasDscal(this->mHBlas, cols,&zero,out+r+start*oSize,oSize), __LINE__);
      }

      GpuMatrix& work = this->cuWorkBuffer();
      work.reshape(r, cols);
      double *pO = out+start*oSize;
      int ldo = oSize;
      CheckCuda(cublasDgeam(this->mHBlas, CUBLAS_OP_N, CUBLAS_OP_N, r, cols, &one, pO, ldo, &zero, pO, ldo, work.data(), r), __LINE__);

      // P to CSR
      this->mcsrP.reshape(r, r, P.size()-1);
      Gpu::makeBidiagCsr<<<1,128,0,this->stream(streamId)>>>(P.data(), r, this->mcsrP.rowPtr(), this->mcsrP.colIdx(), this->mcsrP.data());

      cusparseSpMatDescr_t descrP;
      CheckCuda(cusparseCreateCsr(&descrP, r, r, this->mcsrP.nnz(), this->mcsrP.rowPtr(), this->mcsrP.colIdx(), this->mcsrP.data(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F), __LINE__);

      cusparseDnMatDescr_t descrO;
      CheckCuda(cusparseCreateDnMat(&descrO, r, cols, oSize, out+start*oSize, CUDA_R_64F, CUSPARSE_ORDER_COL), __LINE__);

      cusparseDnMatDescr_t descrT;
      CheckCuda(cusparseCreateDnMat(&descrT, r, cols, r, work.data(), CUDA_R_64F, CUSPARSE_ORDER_COL), __LINE__);

      size_t lworkMM = 0;
      CheckCuda(cusparseSpMM_bufferSize(
            this->mHSparse,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            &one,
            descrP, descrT, &zero, descrO, CUDA_R_64F, CUSPARSE_MM_ALG_DEFAULT, &lworkMM), __LINE__);
      this->mcuBuffer.resize(lworkMM/sizeof(double), GpuMatrix::KEEP_SIZE);

      CheckCuda(cusparseSpMM(
         this->mHSparse,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         &one,
         descrP, 
         descrT,
         &zero,
         descrO, 
         CUDA_R_64F, 
         CUSPARSE_MM_ALG_DEFAULT,
         this->mcuBuffer.data()), __LINE__);

      CheckCuda(cusparseDestroySpMat(descrP), __LINE__);
      CheckCuda(cusparseDestroyDnMat(descrT), __LINE__);
      CheckCuda(cusparseDestroyDnMat(descrO), __LINE__);

      // Set default stream
      this->useStream(-1);
   }
#endif //QUICC_TRIPRODUCT_SPMM

#ifdef QUICC_TRISOLVE_TBSV
   void IWorlandBackend::computeTriSolve(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& S) const
   {
      // Set stream
      this->useStream(streamId);

      int s = S.cols();

      cublasFillMode_t uplo;
      this->getMatrixFillMode(streamId, S, uplo);

      cublasOperation_t trans = CUBLAS_OP_N;
      cublasDiagType_t diag = CUBLAS_DIAG_NON_UNIT;
      int k = 1;
      int lda = 2;
      int incx = 1;

      for(int j = 0; j < cols; j++)
      {
         double *X = out + (start+j)*oSize;
         CheckCuda(cublasDtbsv(this->mHBlas, uplo, trans, diag, s, k, S.data(), lda, X, incx), __LINE__);
      }
   }
#endif //QUICC_TRISOLVE_TBSV

#ifdef QUICC_TRISOLVE_CSRSV2
   void IWorlandBackend::computeTriSolve(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& S) const
   {
      // Set stream
      this->useStream(streamId);

      int s = S.cols();
      const double one = 1.0;

      cublasFillMode_t S_uplo;
      this->getMatrixFillMode(streamId, S, S_uplo);

      cusparseMatDescr_t descrS = nullptr;
      csrsv2Info_t info = nullptr;
      int lworkSV = 0;
      GpuMatrix& work = this->cuWorkBuffer();
      work.reshape(s, 1);

      CheckCuda(cusparseCreateCsrsv2Info(&info), __LINE__);
      CheckCuda(cusparseCreateMatDescr(&descrS), __LINE__);
      CheckCuda(cusparseSetMatIndexBase(descrS,CUSPARSE_INDEX_BASE_ZERO), __LINE__);
      if(S_uplo == CUBLAS_FILL_MODE_LOWER)
      {
         CheckCuda(cusparseSetMatFillMode(descrS, CUSPARSE_FILL_MODE_LOWER), __LINE__);
      } else
      {
         CheckCuda(cusparseSetMatFillMode(descrS, CUSPARSE_FILL_MODE_UPPER), __LINE__);
      }
      CheckCuda(cusparseSetMatDiagType(descrS, CUSPARSE_DIAG_TYPE_NON_UNIT), __LINE__);
      const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_NO_LEVEL;

      this->mcsrS.reshape(s, s, S.size()-1);
      Gpu::makeBidiagCsr<<<1,128,0,this->stream(streamId)>>>(S.data(), s, this->mcsrS.rowPtr(), this->mcsrS.colIdx(), this->mcsrS.data());

      CheckCuda(cusparseDcsrsv2_bufferSize(
         this->mHSparse,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         s,
         this->mcsrS.nnz(),
         descrS,
         this->mcsrS.data(), this->mcsrS.rowPtr(), this->mcsrS.colIdx(),
         info,
         &lworkSV), __LINE__);
      this->mcuBuffer.resize(lworkSV/sizeof(double), GpuMatrix::KEEP_SIZE);

      CheckCuda(cusparseDcsrsv2_analysis(
         this->mHSparse,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         s,
         this->mcsrS.nnz(),
         descrS,
         this->mcsrS.data(), this->mcsrS.rowPtr(), this->mcsrS.colIdx(),
         info,
         policy,
         this->mcuBuffer.data()), __LINE__);

      for(int j = 0; j < cols; j++)
      {
         double *X = out + (start+j)*oSize;
         CheckCuda(cublasDcopyAsync(this->mHBlas, s, X, 1, work.data(), 1, this->stream(streamId)), __LINE__);
         CheckCuda(cusparseDcsrsv2_solve(
            this->mHSparse,
            CUSPARSE_OPERATION_NON_TRANSPOSE,
            s,
            this->mcsrS.nnz(),
            &one,
            descrS,
            this->mcsrS.data(), this->mcsrS.rowPtr(), this->mcsrS.colIdx(),
            info,
            work.data(),
            X,
            policy,
            this->mcuBuffer.data()), __LINE__);
      }

      CheckCuda(cusparseDestroyCsrsv2Info(info), __LINE__);
   }
#endif //QUICC_TRISOLVE_CSRSV2

#ifdef QUICC_TRISOLVE_CSRSM2
   void IWorlandBackend::computeTriSolve(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& S) const
   {
      // Set stream
      this->useStream(streamId);

      int s = S.cols();
      const double one = 1.0;
      const int algo = 1;

      cusparseMatDescr_t descrS = nullptr;
      csrsm2Info_t info = nullptr;
      cublasFillMode_t S_uplo;
      this->getMatrixFillMode(streamId, S, S_uplo);

      CheckCuda(cusparseCreateCsrsm2Info(&info), __LINE__);
      CheckCuda(cusparseCreateMatDescr(&descrS), __LINE__);
      CheckCuda(cusparseSetMatIndexBase(descrS,CUSPARSE_INDEX_BASE_ZERO), __LINE__);
      if(S_uplo == CUBLAS_FILL_MODE_LOWER)
      {
         CheckCuda(cusparseSetMatFillMode(descrS, CUSPARSE_FILL_MODE_LOWER), __LINE__);
      } else
      {
         CheckCuda(cusparseSetMatFillMode(descrS, CUSPARSE_FILL_MODE_UPPER), __LINE__);
      }
      CheckCuda(cusparseSetMatDiagType(descrS, CUSPARSE_DIAG_TYPE_NON_UNIT), __LINE__);
      const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_NO_LEVEL;

      // S to CSR
      this->mcsrS.reshape(s, s, S.size()-1);
      Gpu::makeBidiagCsr<<<1,128,0,this->stream(streamId)>>>(S.data(), s, this->mcsrS.rowPtr(), this->mcsrS.colIdx(), this->mcsrS.data());

      size_t lworkSM = 0;
      CheckCuda(cusparseDcsrsm2_bufferSizeExt(
         this->mHSparse,
         algo,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         s,
         cols,
         this->mcsrS.nnz(),
         &one,
         descrS,
         this->mcsrS.data(), this->mcsrS.rowPtr(), this->mcsrS.colIdx(),
         out+start*oSize,
         oSize,
         info,
         policy,
         &lworkSM), __LINE__);
      this->mcuBuffer.resize(lworkSM/sizeof(double), GpuMatrix::KEEP_SIZE);

      CheckCuda(cusparseDcsrsm2_analysis(
         this->mHSparse,
         algo,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         s,
         cols,
         this->mcsrS.nnz(),
         &one,
         descrS,
         this->mcsrS.data(), this->mcsrS.rowPtr(), this->mcsrS.colIdx(),
         out+start*oSize,
         oSize,
         info,
         policy,
         this->mcuBuffer.data()), __LINE__);

      CheckCuda(cusparseDcsrsm2_solve(
         this->mHSparse,
         algo,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         CUSPARSE_OPERATION_NON_TRANSPOSE,
         s,
         cols,
         this->mcsrS.nnz(),
         &one,
         descrS,
         this->mcsrS.data(), this->mcsrS.rowPtr(), this->mcsrS.colIdx(),
         out+start*oSize,
         oSize,
         info,
         policy,
         this->mcuBuffer.data()), __LINE__);

      CheckCuda(cusparseDestroyCsrsm2Info(info), __LINE__);
      CheckCuda(cusparseDestroyMatDescr(descrS), __LINE__);

      // Set default stream
      this->useStream(-1);
   }
#endif //QUICC_TRISOLVE_CSRSM2

namespace Gpu
{

   __global__ void jacobiShiftMatrix(double* J, const int jSize, const double alpha, const double beta)
   {
      double a1 = alpha + 1.0;
      double ab = alpha + beta;
      double ab1 = ab + 1.0;
      double ab2 = ab + 2.0;
      for(int i = threadIdx.x; i < jSize; i+=blockDim.x)
      {
         double n = static_cast<double>(i);
         J[i] = sqrt(2.0*(n + beta)*(n + ab));
         J[jSize+i] = sqrt((2.0*n + ab)*(2.0*n + ab1));
         J[2*jSize+i] = sqrt(2.0*(n + 1.0)*(n + a1));
         J[3*jSize+i] = sqrt((2.0*n + ab1)*(2.0*n + ab2));
      }

      if(threadIdx.x == 0 && ab == 0.0)
      {
         J[0] = std::sqrt(beta);
         J[jSize] = 1.0;
      }
   }

   __global__ void buildShiftU(double* U, const int i, const int uSize, const double* J, const int jSize, const double norm)
   {
      // Safety asserts
      assert(uSize-1 <= jSize);
      assert(i+uSize-1 <= jSize);
      assert(2*i+uSize <= jSize);
      assert(i+uSize <= jSize);

      for(int j = threadIdx.x; j < uSize-1; j+=blockDim.x)
      {
         U[2*(j+1)] = -J[2*jSize+j]/J[3*jSize+i+j];
      }
      for(int j = threadIdx.x; j < uSize; j+=blockDim.x)
      {
         U[2*j+1] = J[2*i+j]/J[jSize+i+j];
      }

      // Normalize
      if(norm != 1.0)
      {
         for(int j = threadIdx.x; j < uSize-1; j+=blockDim.x)
         {
            U[2*(j+1)] *= norm;
         }
         for(int j = threadIdx.x; j < uSize; j+=blockDim.x)
         {
            U[2*j+1] *= norm;
         }
      }

      // Set upper triangular flag value
      if(threadIdx.x == 0)
      {
         U[0] = IWorlandBackend::UPPER_BANDED;
      }
   }

   __global__ void buildShiftV(double* V, const int i, const int vSize, const double* J, const int jSize, const double norm)
   {
      // Safety asserts
      assert(vSize-1 <= jSize);
      assert(i+vSize-1 <= jSize);
      assert(2*i+vSize <= jSize);
      assert(i+vSize <= jSize);

      for(int j = threadIdx.x; j < vSize-1; j+=blockDim.x)
      {
         V[2*(j+1)] = J[2*jSize+j]/J[3*jSize+i+j];
      }
      for(int j = threadIdx.x; j < vSize; j+=blockDim.x)
      {
         V[2*j+1] = J[2*i+j]/J[jSize+i+j];
      }

      // Normalize
      if(norm != 1.0)
      {
         for(int j = threadIdx.x; j < vSize-1; j+=blockDim.x)
         {
            V[2*(j+1)] *= norm;
         }
         for(int j = threadIdx.x; j < vSize; j+=blockDim.x)
         {
            V[2*j+1] *= norm;
         }
      }

      // Set upper triangular flag value
      if(threadIdx.x == 0)
      {
         V[0]= IWorlandBackend::UPPER_BANDED;
      }
   }

   __global__ void buildShiftM(double* M, const int i, const int mSize, const double* J, const int jSize, const double norm, const bool isSquare)
   {
      // Safety asserts
      assert(2*i+1+mSize <= jSize);
      assert(i+mSize <= jSize);
      assert(mSize-1 <= jSize);
      assert(i+1+mSize-1 <= jSize);

      for(int j = threadIdx.x; j < mSize; j+=blockDim.x)
      {
         M[2*j] = J[2*i+1+j]/J[3*jSize+i+j];
      }
      for(int j = threadIdx.x; j < mSize-1; j+=blockDim.x)
      {
         M[2*j+1] = J[2*jSize+j]/J[jSize+i+1+j];
      }

      // Normalize
      if(norm != 1.0)
      {
         for(int j = threadIdx.x; j < mSize; j+=blockDim.x)
         {
            M[2*j] *= norm;
         }
         for(int j = threadIdx.x; j < mSize-1; j+=blockDim.x)
         {
            M[2*j+1] *= norm;
         }
      }

      // Set lower triangular flag value
      if(threadIdx.x == 0)
      {
         if(isSquare)
         {
            M[2*mSize-1] = IWorlandBackend::LOWER_BANDED;
         } else
         {
            M[2*mSize-1] = IWorlandBackend::LOWER_BANDED_PADDED;
         }
      }
   }

   __global__ void scaleBlock(double* D, const int dSize, const double c, const int rows, const int cols)
   {
      for(int j = 0; j < cols; j++)
      {
         for(int i = threadIdx.x; i < rows; i+=blockDim.x)
         {
            D[i+j*dSize] *= c;
         }
      }
   }

   __global__ void makeBidiagCsr(double* A, const int aSize, int* rowPtr, int* colIdx,  double* val)
   {
      int nnz = 2*aSize - 1;

      if(A[0] == IWorlandBackend::UPPER_BANDED)
      {
         for(int i = threadIdx.x; i < aSize; i+=blockDim.x)
         {
            if(i == aSize-1)
            {
               rowPtr[i] = 2*i;
               colIdx[2*i] = i;
               rowPtr[aSize] = nnz;
            } else
            {
               rowPtr[i] = 2*i;
               colIdx[2*i] = i;
               colIdx[2*i+1] = i+1;
            }
         }

         for(int i = threadIdx.x; i < nnz; i+=blockDim.x)
         {
            val[i] = A[i+1];
         }
      } else if(std::abs(A[2*aSize-1]) == IWorlandBackend::LOWER_BANDED)
      {
         for(int i = threadIdx.x; i < aSize; i+=blockDim.x)
         {
            if(i == 0)
            {
               rowPtr[i] = 0;
               colIdx[2*i] = 0;
               rowPtr[aSize] = nnz;
            } else
            {
               rowPtr[i] = 2*(i-1)+1;
               colIdx[2*i-1] = i-1;
               colIdx[2*i] = i;
            }
         }

         for(int i = threadIdx.x; i < nnz; i+=blockDim.x)
         {
            val[i] = A[i];
         }
      } else
      {
         printf("makeBidiagCsr could not identify if upper or lower matrix should be constructed!");
         asm("trap;");
      }
   }

}

}
}
}
}
}
