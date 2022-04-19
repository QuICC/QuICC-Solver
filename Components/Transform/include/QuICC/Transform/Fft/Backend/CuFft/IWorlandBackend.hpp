/**
 * @file IWorlandBackend.hpp
 * @brief Interface for a generic Worland cuFFT based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_CUFFT_IWORLANDBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_CUFFT_IWORLANDBACKEND_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <set>

// External includes
//
#include "cublas_v2.h"
#include "cusparse.h"

// Project includes
//
#include "QuICC/Transform/Fft/Backend/CuFft/ICuFftBackend.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuMatrix.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuCsrMatrix.hpp"
#include "QuICC/Transform/Fft/Backend/CuFft/GpuPinnedMatrix.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace CuFft {

   /**
    * @brief Interface for a generic Worland cuFFT based integrator
    */
   class IWorlandBackend: public ICuFftBackend
   {
      public:
         /// Typedef for the configuration class
         typedef Worland::Setup SetupType;

         /// Typedef for the configuration class as a shared pointer
         typedef Worland::SharedSetup SharedSetupType;

         /**
          * @brief Typedef for block location
          * <0>: l
          * <1>: start
          * <2>: multiplicity
          * <3>: Accurate resolution
          * <4>: current l
          */
         typedef std::tuple<int,int,int,int,int> LocationType;

         //// Typedef for vector of block locations
         typedef std::vector<LocationType> LocationVector;

         /**
          * @brief Flag for identifying upper triangular matrix
          */
         static const double UPPER_BANDED;

         /**
          * @brief Flag for identifying lower triangular matrix
          */
         static const double LOWER_BANDED;

         /**
          * @brief Flag for identifying lower rectangular triangular matrix
          */
         static const double LOWER_BANDED_PADDED;

         /**
          * @brief Constructor
          */
         IWorlandBackend();

         /**
          * @brief Constructor
          */
         IWorlandBackend(const int nStreams, const int nBatches);

         /**
          * @brief Destructor
          */
         virtual ~IWorlandBackend();

         /**
          * @brief Set zero filter
          */
         void setZFilter(const std::set<int>& filter) const;

         /**
          * @brief Initialise the FFT transforms
          */
         virtual void init(const SetupType& setup, const int lshift, const bool lshiftOnlyParity = false, const bool alwaysZeroNegative = false) const;

         /**
          * @brief Change internal spectral resolution
          */
         virtual void setWSize(const unsigned int shiftMaxL) const;

         /**
          * @brief Add temporary storage
          */
         void addStorage(const int inExtras, const int outExtras) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void io(double* out, const double* in) const;

         /**
          * @brief Scaling by constant
          */
         void scaleC(const double c, const bool isEven, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by a*l + y
          */
         void scaleALPY(const double a, const double y, const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by 2.0(l+n+1)*sqrt((n+1)/(l+n+1))
          */
         void scaleD(const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by ??
          */
         void scaleSphLaplA(const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by ??
          */
         void scaleSphLaplB(const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

         /**
          * @brief Shift l in temporary storage
          */
         void lshift(const unsigned int id, const int lshift, const bool isEven) const;

         /**
          * @brief Shift n expansion in temporary storage
          */
         void nshift(const unsigned int id, const int nshift, const bool isEven) const;

         /**
          * @brief Copy out temporary storage into new extra temporary
          */
         void copy(const int to, const int from, const int nshift, const bool isEven) const;

         /**
          * @brief Add extra temporary to out temporary storage
          */
         void add(const int to, const int from, const int nshift, const bool isEven) const;

      protected:
         /**
          * @brief Is physical representation even?
          */
         bool isPhysEven(const bool isSpecEven) const;

         /**
          * @brief Is filtered out by zero filter
          */
         bool inZFilter(const int l) const;

         /**
          * @brief Get blocksize
          */
         int blockSize(const bool isEven) const;

         /**
          * @brief Accurate expansion for given l
          */
         virtual int lSize(const int l) const = 0;

         /**
          * @brief Accurate expansion for given index i and l
          */
         virtual int lSize(const int i, const int l) const = 0;

         /**
          * @brief Get pointer to vector of locations
          */
         LocationVector* pLoc(const bool isEven, const unsigned int id = 0) const;

         /**
          * @brief Reset locations to initial values
          */
         void resetLocations(const bool isEven, const int id = 0) const;

         /**
          * @brief Initialize temporary work storage
          */
         void initStorage(const int rows, const int cols, const int n, std::vector<GpuMatrix>& s) const;

         /**
          * @brief Set even or odd plan
          */
         void setPlan(const bool isEven) const;

         /**
          * @brief Initialize J matrices
          */
         void initJ() const;

         /**
          * @brief Get J matrix and size based on parity flag
          */
         GpuMatrix& J(const bool isEven) const;

         /**
          * @brief Get work temporary
          */
         GpuMatrix& cuWorkTmp(const unsigned int id = 0) const;

         /**
          * @brief Get work buffer temporary
          */
         GpuMatrix& cuWorkBuffer() const;


         /**
          * @brief Build U(upper) for shifting jacobi alpha by +1
          */
         void buildShiftU(GpuMatrix& U, const double alpha, const double beta, const double norm = 1.0) const;


         /**
          * @brief Build V(upper) for shifting jacobi beta by +1
          */
         void buildShiftV(GpuMatrix& V, const double alpha, const double beta, const double norm = 1.0) const;


         /**
          * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
          */
         void buildShiftM(GpuMatrix& M, const double alpha, const double beta, const double norm = 1.0, const bool isSquare = true) const;


         /**
          * @brief Build V(upper),M(lower) pair for shifting jacobi beta by +1
          */
         void buildShiftPair(const int streamId, GpuMatrix& V, GpuMatrix& M, const int i, const GpuMatrix& J, const double normV = 1.0, const double normM = 1.0, const bool isSquare = true) const;

         /**
          * @brief Apply triangular banded product
          */
         void applyTriProduct(double* out, const int oSize, const int start, const int cols, GpuMatrix& A) const;

         /**
          * @brief Apply triangular banded product
          */
         void applyBandProduct(double* out, const int oSize, const int start, const int cols, GpuMatrix& A) const;

         /**
          * @brief Apply triangular banded linear solve
          */
         void applyTriSolve(double* out, const int oSize, const int start, const int cols, GpuMatrix& A) const;

         /**
          * @brief Apply pair of banded product and banded linear solve
          */
         void applyPair(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& P, GpuMatrix& S) const;

         /**
          * @brief Get triangular banded storage fill type
          */
         bool getMatrixFillMode(const int streamId, const GpuMatrix& A, cublasFillMode_t& uplo) const;

         /**
          * @brief Get triangular banded storage fill type
          */
         void getMatrixBands(const GpuMatrix& A, int& kl, int& ku) const;

         /**
          * @brief Compute triangular banded product
          */
         void computeTriProduct(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& A) const;

         /**
          * @brief Compute triangular banded linear solve
          */
         void computeTriSolve(const int streamId, double* out, const int oSize, const int start, const int cols, GpuMatrix& A) const;

         /**
          * @brief Use given stream for CUBLAS and CUSPARSE
          */
         void useStream(const int id) const;

         /**
          * @brief Use given stream
          */
         void useStream(cublasHandle_t h, const int id) const;

         /**
          * @brief Use given stream
          */
         void useStream(cusparseHandle_t h, const int id) const;

         /**
          * @brief Minimum number of blocks per batch
          */
         static const int MIN_BATCH_BLOCKSIZE;

         /**
          * @brief CUBLAS handle
          */
         mutable cublasHandle_t mHBlas;

         /**
          * @brief CUSPARSE handle
          */
         mutable cusparseHandle_t mHSparse;

         /**
          * @brief Temporary GPU forward data
          */
         mutable std::vector<cufftDoubleReal*>  mcuFwd;

         /**
          * @brief Temporary GPU backward data
          */
         mutable std::vector<cufftDoubleComplex*>  mcuBwd;

         /**
          * @brief Temporary GPU real work data
          */
         mutable std::vector<cufftDoubleReal*>  mcuWork;

         /**
          * @brief Temporary GPU complex work data
          */
         mutable std::vector<cufftDoubleComplex*>  mcuWorkZ;

         /**
          * @brief Even harmonic degree block size
          */
         mutable std::vector<int> mcuEBlockSize;

         /**
          * @brief Odd harmonic degree block size
          */
         mutable std::vector<int> mcuOBlockSize;

         /**
          * @brief Plan switch
          */
         mutable bool mPlanIsEven;

         /**
          * @brief Flip parity with respect to l
          */
         mutable bool mFlipped;

         /**
          * @brief Spec size
          */
         mutable int mSpecSize;

         /**
          * @brief Forward FFT size
          */
         mutable int mFwdSize;

         /**
          * @brief Backward FFT size
          */
         mutable int mBwdSize;

         /**
          * @brief Biggest Worland expansion size (ie. l = 0)
          */
         mutable int mWSize;

         /**
          * @brief Even harmonic degree block size
          */
         mutable int mEBlockSize;

         /**
          * @brief Odd harmonic degree block size
          */
         mutable int mOBlockSize;

         /**
          * @brief Size of J for even modes
          */
         mutable int mJEvenSize;

         /**
          * @brief Size of J for odd modes
          */
         mutable int mJOddSize;

         /**
          * @brief Filter to zero modes
          */
         mutable std::set<int> mZFilter;

         /**
          * @brief Location of even harmonic degrees
          */
         mutable std::vector<LocationVector> mELoc;

         /**
          * @brief Location of odd harmonic degrees
          */
         mutable std::vector<LocationVector> mOLoc;

         /**
          * @brief Location of zero-ed harmonic degrees
          */
         mutable LocationVector mZLoc;

         /**
          * @brief Jacobi shift matrix even modes
          */
         mutable GpuMatrix mcuJEven;

         /**
          * @brief Jacobi shift matrix odd modes
          */
         mutable GpuMatrix mcuJOdd;

         /**
          * @brief GPU storage for U matrix
          */
         mutable GpuMatrix mcuU;

         mutable GpuCsrMatrix mcsrP;
         mutable GpuCsrMatrix mcsrS;
         mutable GpuMatrix mcuBuffer;

         /**
          * @brief GPU storage for V matrix
          */
         mutable GpuMatrix mcuV;

         /**
          * @brief GPU storage for M matrix
          */
         mutable GpuMatrix mcuM;

         /**
          * @brief Host pinned staging storage
          */
         mutable GpuPinnedMatrix mPinnedIn;

         /**
          * @brief Host pinned staging storage
          */
         mutable GpuPinnedMatrix mPinnedOut;

         /**
          * @brief Temporary GPU data
          */
         mutable std::vector<GpuMatrix>  mcuInTmp;

         /**
          * @brief Temporary GPU data
          */
         mutable std::vector<GpuMatrix>  mcuOutTmp;

         /**
          * @brief Pointer to temporary data
          */
         mutable std::vector<GpuMatrix>*  mpcuWorkTmp;

         /**
          * @brief Input data pointer
          */
         mutable const double* mpIn;

         /**
          * @brief Out data pointer
          */
         mutable double* mpOut;

      private:
   };

   namespace Gpu
   {
      /**
       * @brief Compute shift matrix for Jacobi polynomials: P^{a,b-1} -> P^{a,b}
       */
      __global__ void jacobiShiftMatrix(double* J, const int jSize, const double alpha, const double beta);

      /**
       * @brief Build U(upper) for shifting jacobi alpha by +1
       */
      __global__ void buildShiftU(double* U, const int i, const int uSize, const double* J, const int jSize, const double norm = 1.0);

      /**
       * @brief Build V(upper) for shifting jacobi beta by +1
       */
      __global__ void buildShiftV(double* V, const int i, const int vSize, const double* J, const int jSize, const double norm = 1.0);

      /**
       * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
       */
      __global__ void buildShiftM(double* M, const int i, const int mSize, const double* J, const int jSize, const double norm = 1.0, const bool isSquare = true);

      /**
       * @brief Scale block data
       */
      __global__ void scaleBlock(double* D, const int dSize, const double c, const int rows, const int cols);

      /**
       * @brief Convert BLAS banded storage to CSR
       */
      __global__ void makeBidiagCsr(double* A, const int aSize, int* rowPtr, int* colInd, double* val);
   }

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_CUFFT_IWORLANDBACKEND_HPP
