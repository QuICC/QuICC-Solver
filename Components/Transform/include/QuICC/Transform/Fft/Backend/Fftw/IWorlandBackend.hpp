/**
 * @file IWorlandBackend.hpp
 * @brief Interface for a generic Worland FFTW based integrator
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_IWORLANDBACKEND_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_IWORLANDBACKEND_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Transform/Fft/Backend/Fftw/IFftwBackend.hpp"
#include "QuICC/Transform/Fft/Worland/Setup.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Interface for a generic Worland FFTW based integrator
    */
   class IWorlandBackend: public IFftwBackend
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
          * @brief Constructor
          */
         IWorlandBackend();

         /**
          * @brief Destructor
          */
         virtual ~IWorlandBackend();

         /**
          * @brief Set zero filter
          */
         void setZFilter(const std::set<int>& filter) const;

         /**
          * @brief Initialise the FFTW transforms
          */
         virtual void init(const SetupType& setup, const int lshift, const int extraN, const bool lshiftOnlyParity = false, const bool alwaysZeroNegative = false) const;

         /**
          * @brief Initialise the FFTW transforms
          */
         void addStorage(const int inExtras, const int outExtras) const;

         /**
          * @brief Set input and output data pointers for FFT (R2R)
          */
         virtual void io(MHDFloat* out, const MHDFloat* in) const;

         /**
          * @brief Scaling by constant
          */
         void scaleC(const MHDFloat c, const bool isEven, const unsigned int id = 0) const;

         /**
          * @brief Scale coefficients by a*l + y
          */
         void scaleALPY(const MHDFloat a, const MHDFloat y, const bool isEven, const int lshift = 0, const unsigned int id = 0) const;

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
          * @brief Change internal spectral resolution
          */
         virtual void setWSize() const;

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
         void initStorage(const int rows, const int cols, const int n, std::vector<Matrix>& s) const;

         /**
          * @brief Set even or odd plan
          */
         void setPlan(const bool isEven) const;

         /**
          * @brief Initialize J matrices
          */
         void initJ() const;

         /**
          * @brief Get J matrix based on parity flag
          */
         const Matrix& J(const bool isEven) const;

         /**
          * @brief Get banded matrices based on parity flag
          */
         const std::vector<Matrix>& banded(const bool isEven) const;

         /**
          * @brief Get work temporary
          */
         Matrix& workTmp(const unsigned int id = 0) const;

         /**
          * @brief Convert banded storage upper matrix to sparse matrix
          */
         void makeUpper(SparseMatrix& M, const Matrix& Mp) const;

         /**
          * @brief Convert banded storage lower matrix to sparse matrix
          */
         void makeLower(SparseMatrix& M, const Matrix& Mp, const bool isSquare = true) const;

         /**
          * @brief Compute shift matrix for Jacobi polynomials: P^{a,b-1} -> P^{a,b}
          */
         void jacobiShiftMatrix(Matrix& rJ, const int maxN, const MHDFloat alpha, const MHDFloat beta) const;


         /**
          * @brief Build U(upper) for shifting jacobi alpha by +1
          */
         void buildShiftU(Matrix& U, const int maxN, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm = 1.0) const;

         /**
          * @brief Build U(upper) for shifting jacobi alpha by +1
          */
         void buildShiftU(SparseMatrix& U, const int maxN, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm = 1.0) const;

         /**
          * @brief Build U(upper) for shifting jacobi alpha by +1
          */
         void buildShiftU(Matrix& U, const int i, const int maxN, const Matrix& J, const MHDFloat norm = 1.0) const;

         /**
          * @brief Build U(upper) for shifting jacobi alpha by +1
          */
         void buildShiftU(SparseMatrix& U, const int i, const int maxN, const Matrix& J, const MHDFloat norm = 1.0) const;


         /**
          * @brief Build V(upper) for shifting jacobi beta by +1
          */
         void buildShiftV(Matrix& V, const int maxN, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm = 1.0) const;

         /**
          * @brief Build V(upper) for shifting jacobi beta by +1
          */
         void buildShiftV(SparseMatrix& V, const int maxN, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm = 1.0) const;

         /**
          * @brief Build V(upper) for shifting jacobi beta by +1
          */
         void buildShiftV(Matrix& V, const int i, const int maxN, const Matrix& J, const MHDFloat norm = 1.0) const;

         /**
          * @brief Build V(upper) for shifting jacobi beta by +1
          */
         void buildShiftV(SparseMatrix& V, const int i, const int maxN, const Matrix& J, const MHDFloat norm = 1.0) const;


         /**
          * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
          */
         void buildShiftM(Matrix& M, const int maxN, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm = 1.0, const bool isSquare = true) const;

         /**
          * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
          */
         void buildShiftM(SparseMatrix& M, const int maxN, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm = 1.0, const bool isSquare = true) const;

         /**
          * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
          */
         void buildShiftM(Matrix& M, const int i, const int maxN, const Matrix& J, const MHDFloat norm = 1.0, const bool isSquare = true) const;

         /**
          * @brief Build M(lower) for shifting jacobi beta by +1 and multiply by (1+x)
          */
         void buildShiftM(SparseMatrix& M, const int i, const int maxN, const Matrix& J, const MHDFloat norm = 1.0, const bool isSquare = true) const;


         /**
          * @brief Build V(upper),M(lower) pair for shifting jacobi beta by +1
          */
         void buildShiftPair(Matrix& PS, const bool isVMOrder, const int i, const int maxN, const Matrix& J, const MHDFloat normV = 1.0, const MHDFloat normM = 1.0, const bool isSquare = true) const;

         /**
          * @brief Build V(upper),M(lower) pair for shifting jacobi beta by +1
          */
         void buildShiftPair(SparseMatrix& V, SparseMatrix& M, const int i, const int maxN, const Matrix& J, const MHDFloat normV = 1.0, const MHDFloat normM = 1.0, const bool isSquare = true) const;

         /**
          * @brief Apply triangular banded product
          */
         void applyTriProduct(Matrix& out, const int start, const int cols, Matrix& A) const;

         /**
          * @brief Apply triangular sparse product
          */
         void applyTriProduct(Matrix& out, const int start, const int cols, const SparseMatrix& A) const;

         /**
          * @brief Apply triangular banded product
          */
         void applyBandProduct(Matrix& out, const int start, const int cols, Matrix& A) const;

         /**
          * @brief Apply triangular sparse product
          */
         void applyBandProduct(Matrix& out, const int start, const int cols, const SparseMatrix& A) const;

         /**
          * @brief Apply triangular banded linear solve
          */
         void applyTriSolve(Matrix& out, const int start, const int cols, Matrix& A) const;

         /**
          * @brief Apply triangular sparse linear solve
          */
         void applyTriSolve(Matrix& out, const int start, const int cols, const SparseMatrix& A) const;

         /**
          * @brief Apply pair of banded product and banded linear solve
          */
         void applyPair(Matrix& out, const int start, const int cols, const int rows, const Matrix& PS) const;

         /**
          * @brief Apply pair of sparse product and sparse linear solve
          */
         void applyPair(Matrix& out, const int start, const int cols, const int rows, const SparseMatrix& P, const SparseMatrix& S) const;

         /**
          * @brief Flag for identifying upper triangular matrix
          */
         static const MHDFloat UPPER_BANDED;

         /**
          * @brief Flag for identifying lower triangular matrix
          */
         static const MHDFloat LOWER_BANDED;

         /**
          * @brief Flag for identifying lower rectangular triangular matrix
          */
         static const MHDFloat LOWER_BANDED_PADDED;

         /**
          * @brief Flip parity with respect to l
          */
         mutable bool mFlipped;

         /**
          * @brief Plan for the transform of odd harmonic degree
          */
         mutable fftw_plan*   mpPlan;

         /**
          * @brief Plan for the transform of odd harmonic degree
          */
         mutable fftw_plan   mOddPlan;

         /**
          * @brief Spec size
          */
         mutable int mSpecSize;

         /**
          * @brief For some operators, the truncation needs to be extended by a few extra modes i
          *        in order to produce mSpecSize accurate modes at the end. For example, the use
          *        of the second order quas-inverse I2 requires 3 additional modes due to the 3
          *        super diagonals.
          */
         mutable int mWExtra;

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
         mutable Matrix mJEven;

         /**
          * @brief Jacobi shift matrix odd modes
          */
         mutable Matrix mJOdd;

         /**
          * @brief Odd l Banded matrices
          */
         mutable std::vector<Matrix> mOBanded;

         /**
          * @brief Even l Banded matrices
          */
         mutable std::vector<Matrix> mEBanded;

         /**
          * @brief Pointer to temporary data
          */
         mutable std::vector<Matrix>*  mpWorkTmp;

         /**
          * @brief Temporary data
          */
         mutable std::vector<Matrix>  mInTmp;

         /**
          * @brief Temporary data
          */
         mutable std::vector<Matrix>  mOutTmp;

         /**
          * @brief Input data pointer
          */
         mutable const MHDFloat* mpIn;

         /**
          * @brief Out data pointer
          */
         mutable MHDFloat* mpOut;

      private:
   };

   /**
    * @brief Apply pair of banded upper product and banded lower linear solve
    */
   void applyUpperLower(double* X, const int q0, const int r, const Matrix& PS);

   /**
    * @brief Apply pair of banded lower product and banded upper linear solve
    */
   void applyLowerUpper(double* X, const int q0, const int r, const Matrix& PS);

   /**
    * @brief Apply pair of banded upper product and banded lower linear solve
    */
   void applyUpperLower(Eigen::Ref<Matrix> out, const int r, const Matrix& PS);

   /**
    * @brief Apply pair of banded lower product and banded upper linear solve
    */
   void applyLowerUpper(Eigen::Ref<Matrix> out, const int r, const Matrix& PS);

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_IWORLANDBACKEND_HPP
