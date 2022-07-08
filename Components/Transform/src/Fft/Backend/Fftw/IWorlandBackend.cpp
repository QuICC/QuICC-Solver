/**
 * @file IWorlandBackend.cpp
 * @brief Source of the interface for a generic FFTW based Worland integrator
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/src/misc/blas.h>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/IWorlandBackend.hpp"

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"

#define QUICC_IWORLANDBACKEND_BLAS 0
#define QUICC_IWORLANDBACKEND_SCALAR 1
#define QUICC_IWORLANDBACKEND_VECTOR 2
#if defined QUICC_WORLAND_BACKEND_EIGEN
   #define QUICC_IWORLANDBACKEND QUICC_IWORLANDBACKEND_VECTOR
#elif defined QUICC_WORLAND_BACKEND_BLAS
   #define QUICC_IWORLANDBACKEND QUICC_IWORLANDBACKEND_BLAS
#else
   #define QUICC_IWORLANDBACKEND QUICC_IWORLANDBACKEND_SCALAR
#endif

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   const MHDFloat IWorlandBackend::UPPER_BANDED = 424242.424242;

   const MHDFloat IWorlandBackend::LOWER_BANDED = -424242.424242;

   const MHDFloat IWorlandBackend::LOWER_BANDED_PADDED = -424242.434343;

   IWorlandBackend::IWorlandBackend()
      : mFlipped(false)
   {
      this->mELoc.push_back(LocationVector());
      this->mOLoc.push_back(LocationVector());
   }

   IWorlandBackend::~IWorlandBackend()
   {
   }

   void IWorlandBackend::initStorage(const int rows, const int cols, const int n, std::vector<Matrix>& s) const
   {
      assert(rows > 0);
      assert(cols > 0);
      assert(n > 0);

      // Input temporary storage
      s.reserve(n);
      s.push_back(Matrix(rows, cols));
      s.at(0).setZero(rows, cols);

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
      this->mInTmp.reserve(this->mInTmp.size() + inExtras);
      for(int i = 0; i < inExtras; i++)
      {
         this->mInTmp.push_back(this->mInTmp.at(0));
      }

      // Initialize additional put temporary storage
      this->mOutTmp.reserve(this->mOutTmp.size() + outExtras);
      for(int i = 0; i < outExtras; i++)
      {
         this->mOutTmp.push_back(this->mOutTmp.at(0));
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
      int lt = 0;
      int lf = 0;
      if(this->pLoc(true)->rbegin() != this->pLoc(true)->rend())
      {
         lt = std::get<0>(*this->pLoc(true)->rbegin());
      }
      if(this->pLoc(false)->rbegin() != this->pLoc(false)->rend())
      {
         lf = std::get<0>(*this->pLoc(false)->rbegin());
      }
      // Triangular truncation requirements: WSize = n_lmax_modes + lmax/2
      this->mWSize = this->mSpecSize + this->mWExtra + (std::max(lt,lf))/2;
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

   Matrix& IWorlandBackend::workTmp(const unsigned int id) const
   {
      assert(this->mpWorkTmp);
      assert(id >= 0);
      assert(this->mpWorkTmp->size() > id);

      return this->mpWorkTmp->at(id);
   }

   void IWorlandBackend::initJ() const
   {
      MHDFloat alpha = -0.5;

      // Even modes
      if(this->mELoc.at(0).size() > 0)
      {
         // Max size is L + N + 1
         int lmax = std::get<0>(*this->mELoc.at(0).rbegin());
         int jSize = this->lSize(lmax) + lmax  + 1;
         MHDFloat beta = 0.5;
         this->jacobiShiftMatrix(this->mJEven, jSize, alpha, beta);
      }

      // Odd Modes
      if(this->mOLoc.at(0).size() > 0)
      {
         // Max size is L + N + 1
         int lmax = std::get<0>(*this->mOLoc.at(0).rbegin());
         int jSize = this->lSize(lmax) + lmax  + 1;
         MHDFloat beta = 1.5;
         this->jacobiShiftMatrix(this->mJOdd, jSize, alpha, beta);
      }
   }

   const Matrix& IWorlandBackend::J(const bool isEven) const
   {
      if(this->isPhysEven(isEven))
      {
         return this->mJEven;
      } else
      {
         return this->mJOdd;
      }
   }

   const std::vector<Matrix>& IWorlandBackend::banded(const bool isEven) const
   {
      if(this->isPhysEven(isEven))
      {
         return this->mEBanded;
      } else
      {
         return this->mOBanded;
      }
   }

   void IWorlandBackend::makeUpper(SparseMatrix& M, const Matrix& Mp) const
   {
      int n = Mp.cols();
      M.resize(n, n);
      M.reserve(ArrayI::Constant(n,2));
      M.insert(0,0) = Mp(1,0);
      for(int i = 1; i < n; ++i)
      {
         M.insert(i-1,i) = Mp(0,i);
         M.insert(i,i) = Mp(1,i);
      }
      M.makeCompressed();
   }

   void IWorlandBackend::makeLower(SparseMatrix& M, const Matrix& Mp, const bool isSquare) const
   {
      int n = Mp.cols();
      M.resize(n, n - static_cast<int>(!isSquare));
      M.reserve(ArrayI::Constant(n,2));
      for(int i = 0; i < n-1; ++i)
      {
         M.insert(i,i) = Mp(0,i);
         M.insert(i+1,i) = Mp(1,i);
      }
      if(isSquare)
      {
         M.insert(n-1,n-1) = Mp(0,n-1);
      }
      M.makeCompressed();
   }

   void IWorlandBackend::jacobiShiftMatrix(Matrix& rJ, const int size, const MHDFloat alpha, const MHDFloat beta) const
   {
      MHDFloat a1 = alpha + 1.0;
      MHDFloat ab = alpha + beta;
      MHDFloat ab1 = ab + 1.0;
      MHDFloat ab2 = ab + 2.0;
      Array n = Array::LinSpaced(size, 0.0, static_cast<MHDFloat>(size-1));
      rJ.resize(size,4);
      rJ.col(0) = (2.0*(n.array() + beta)*(n.array() + ab)).sqrt();
      rJ.col(1) = ((2.0*n.array() + ab)*(2.0*n.array() + ab1)).sqrt();
      if(ab == 0.0)
      {
         rJ.col(0)(0) = std::sqrt(beta);
         rJ.col(1)(0) = 1.0;
      }
      rJ.col(2) = (2.0*(n.array() + 1.0)*(n.array() + a1)).sqrt();
      rJ.col(3) = ((2.0*n.array() + ab1)*(2.0*n.array() + ab2)).sqrt();
      if(ab1 == 0.0)
      {
         rJ.col(0)(1) = 1.0;
         rJ.col(3)(0) = 1.0;
      }
   }

   void IWorlandBackend::buildShiftU(Matrix& U, const int size, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm) const
   {
      Matrix J;
      this->jacobiShiftMatrix(J, size, beta, alpha);
      this->buildShiftU(U, 0, size, J, norm);
   }

   void IWorlandBackend::buildShiftU(SparseMatrix& U, const int size, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm) const
   {
      Matrix Up;
      this->buildShiftU(Up, size, alpha, beta, norm);
      this->makeUpper(U, Up);
   }

   void IWorlandBackend::buildShiftU(Matrix& U, const int i, const int size, const Matrix& J, const MHDFloat norm) const
   {
      const int n = size;
      U.resize(2, n);

      // Safety asserts
      assert(n-1 <= J.rows());
      assert(i+n-1 <= J.rows());
      assert(2*i+n <= J.rows());
      assert(i+n <= J.rows());

      Eigen::Map<const Array> Voff_num(J.data()+2*J.rows(), n-1);
      Eigen::Map<const Array> Voff_den(J.data()+3*J.rows()+i, n-1);
      U.row(0).rightCols(n-1) = -Voff_num.array()/Voff_den.array();
      Eigen::Map<const Array> Vdia_num(J.data()+2*i, n);
      Eigen::Map<const Array> Vdia_den(J.data()+J.rows()+i, n);
      U.row(1) = Vdia_num.array()/Vdia_den.array();

      // Normalize
      if(norm != 1.0)
      {
         U *= norm;
      }

      // Set upper triangular flag value
      U(0,0) = UPPER_BANDED;
   }

   void IWorlandBackend::buildShiftU(SparseMatrix& U, const int i, const int size, const Matrix& J, const MHDFloat norm) const
   {
      Matrix Up;
      this->buildShiftU(Up, i, size, J, norm);
      this->makeUpper(U, Up);
   }

   void IWorlandBackend::buildShiftV(Matrix& V, const int size, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm) const
   {
      Matrix J;
      this->jacobiShiftMatrix(J, size, alpha, beta);
      this->buildShiftV(V, 0, size, J, norm);
   }

   void IWorlandBackend::buildShiftV(SparseMatrix& V, const int size, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm) const
   {
      Matrix Vp;
      this->buildShiftV(Vp, size, alpha, beta, norm);
      this->makeUpper(V, Vp);
   }

   void IWorlandBackend::buildShiftV(Matrix& V, const int i, const int size, const Matrix& J, const MHDFloat norm) const
   {
      const int n = size;
      V.resize(2,n);

      // Safety asserts
      assert(n-1 <= J.rows());
      assert(i+n-1 <= J.rows());
      assert(2*i+n <= J.rows());
      assert(i+n <= J.rows());

      Eigen::Map<const Array> Voff_num(J.data()+2*J.rows(), n-1);
      Eigen::Map<const Array> Voff_den(J.data()+3*J.rows()+i, n-1);
      V.row(0).rightCols(n-1) = Voff_num.array()/Voff_den.array();
      Eigen::Map<const Array> Vdia_num(J.data()+2*i, n);
      Eigen::Map<const Array> Vdia_den(J.data()+J.rows()+i, n);
      V.row(1) = Vdia_num.array()/Vdia_den.array();

      // Normalize
      if(norm != 1.0)
      {
         V *= norm;
      }

      // Set upper triangular flag value
      V(0,0) = UPPER_BANDED;
   }

   void IWorlandBackend::buildShiftV(SparseMatrix& V, const int i, const int size, const Matrix& J, const MHDFloat norm) const
   {
      Matrix Vp;
      this->buildShiftV(Vp, i, size, J, norm);
      this->makeUpper(V, Vp);
   }

   void IWorlandBackend::buildShiftM(Matrix& M, const int size, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm, const bool isSquare) const
   {
      Matrix J;
      this->jacobiShiftMatrix(J, size+1, alpha, beta-1);
      this->buildShiftM(M, 0, size, J, norm, isSquare);
   }

   void IWorlandBackend::buildShiftM(SparseMatrix& M, const int size, const MHDFloat alpha, const MHDFloat beta, const MHDFloat norm, const bool isSquare) const
   {
      Matrix Mp;
      this->buildShiftM(Mp, size, alpha, beta, norm, isSquare);
      this->makeLower(M, Mp, isSquare);
   }

   void IWorlandBackend::buildShiftM(Matrix& M, const int i, const int size, const Matrix& J, const MHDFloat norm, const bool isSquare) const
   {
      const int n = size;
      M.resize(2, n);
      M.setConstant(-4242.4242);

      // Safety asserts
      assert(2*i+1+n <= J.rows());
      assert(i+n <= J.rows());
      assert(n-1 <= J.rows());
      assert(i+1+n-1 <= J.rows());

      Eigen::Map<const Array> Mdia_num(J.data()+2*i+1, n);
      Eigen::Map<const Array> Mdia_den(J.data()+3*J.rows()+i, n);
      M.row(0) = Mdia_num.array() / Mdia_den.array();
      Eigen::Map<const Array> Moff_num(J.data()+2*J.rows(), n-1);
      Eigen::Map<const Array> Moff_den(J.data()+1*J.rows()+i+1, n-1);
      M.row(1).leftCols(n-1) = Moff_num.array()/Moff_den.array();

      // Normalize
      if(norm != 1.0)
      {
         M *= norm;
      }

      // Set lower triangular flag value
      if(isSquare)
      {
         M(1,n-1) = LOWER_BANDED;
      } else
      {
         M(1,n-1) = LOWER_BANDED_PADDED;
      }
   }

   void IWorlandBackend::buildShiftM(SparseMatrix& M, const int i, const int size, const Matrix& J, const MHDFloat norm, const bool isSquare) const
   {
      Matrix Mp;
      this->buildShiftM(Mp, i, size, J, norm);
      this->makeLower(M, Mp, isSquare);
   }

   void IWorlandBackend::buildShiftPair(Matrix& PS, const bool isVMOrder, const int i, const int size, const Matrix& J, const MHDFloat normV, const MHDFloat normM, const bool isSquare) const
   {
      Matrix V;
      this->buildShiftV(V, i, size, J, normV);

      Matrix M;
      this->buildShiftM(M, i, size, J, normM, isSquare);

#if QUICC_IWORLANDBACKEND!=QUICC_IWORLANDBACKEND_BLAS
      int cols = V.cols();
      PS.resize(3, cols);
      if(isVMOrder)
      {
         PS.row(0) = V.row(1);
         PS.block(1,0,1, cols-1) = V.block(0,1,1, cols-1);
         PS.block(2,1,1, cols-1) = -M.block(1,0,1, cols-1);
         PS.row(0).array() /= M.row(0).array();
         PS.row(1).array() /= M.row(0).array();
         PS.row(2).array() /= M.row(0).array();
         PS(2,0) = V(0,0);
      } else
      {
         PS.block(0,1,1, cols-1) = M.block(1,0,1, cols-1);
         PS.row(1) = M.row(0);
         PS.block(2,0,1, cols-1) = -V.block(0,1,1, cols-1);
         PS.row(0).array() /= V.row(1).array();
         PS.row(1).array() /= V.row(1).array();
         PS.row(2).array() /= V.row(1).array();
         PS(2,cols-1) = M(1,cols-1);
      }
#else
      PS.resize(4,V.cols());
      if(isVMOrder)
      {
         PS.topRows(2) = V;
         PS.bottomRows(2) = M;
      } else
      {
         PS.topRows(2) = M;
         PS.bottomRows(2) = V;
      }
#endif

   }

   void IWorlandBackend::buildShiftPair(SparseMatrix& V, SparseMatrix& M, const int i, const int size, const Matrix& J, const MHDFloat normV, const MHDFloat normM, const bool isSquare) const
   {
      this->buildShiftV(V, i, size, J, normV);

      this->buildShiftM(M, i, size, J, normM, isSquare);
   }

   void IWorlandBackend::setPlan(const bool isEven) const
   {
      if(this->isPhysEven(isEven))
      {
         this->mpPlan = &this->mPlan;
      } else
      {
         this->mpPlan = &this->mOddPlan;
      }
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

   void IWorlandBackend::applyTriProduct(Matrix& out, const int start, const int cols, const SparseMatrix& A) const
   {
      int r = A.rows();
      int c = A.cols();
      out.block(0, start, r, cols) = A*out.block(0, start, c, cols);
   }

   void IWorlandBackend::applyTriProduct(Matrix& out, const int start, const int cols, Matrix& A) const
   {
      int r = A.cols();
      char UPLO;
      if(A(0,0) == UPPER_BANDED)
      {
         UPLO = 'U';
      } else
      {
         assert(A(1,A.cols()-1) == LOWER_BANDED || A(1,A.cols()-1) == LOWER_BANDED_PADDED);
         UPLO = 'L';
         if(A(1,A.cols()-1) == LOWER_BANDED_PADDED)
         {
            out.block(r, start, 1, cols).setZero();
         }
      }
      char TRANS = 'N';
      char DIAG = 'N';
      int K = 1;
      int LDA = 2;
      int INCX = 1;
      for(int j = 0; j < cols; j++)
      {
         double *X = out.data() + (start+j)*out.rows();
         BLASFUNC(dtbmv)(&UPLO, &TRANS, &DIAG, &r, &K, A.data(), &LDA, X, &INCX);
      }
   }

   void IWorlandBackend::applyBandProduct(Matrix& out, const int start, const int cols, const SparseMatrix& A) const
   {
      int r = A.rows();
      int c = A.cols();
      out.block(0, start, r, cols) = A*out.block(0, start, c, cols);
   }

   void IWorlandBackend::applyBandProduct(Matrix& out, const int start, const int cols, Matrix& A) const
   {
      int r = A.cols();
      Matrix tmp = out.block(0, start, r, cols);
      int KL = A.bottomRightCorner(1,1)(0,0);
      int KU = A.topLeftCorner(1,1)(0,0);

      char TRANS = 'N';
      double ALPHA = 1.0;
      int LDA = KL + KU + 1;
      int INCX = 1;
      double BETA = 0.0;
      for(int j = 0; j < cols; j++)
      {
         double *X = tmp.data() + j*r;
         double *Y = out.data() + (start+j)*out.rows();
         BLASFUNC(dgbmv)(&TRANS, &r, &r, &KL, &KU, &ALPHA, A.data(), &LDA, X, &INCX, &BETA, Y, &INCX);
      }
   }

   void IWorlandBackend::applyTriSolve(Matrix& out, const int start, const int cols, const SparseMatrix& A) const
   {
      int r = A.cols();
      Eigen::SparseLU<SparseMatrix> solver;
      solver.compute(A);
      Matrix tmp = out.block(0, start, r, cols);
      Matrix tmp2 = solver.solve(tmp);
      out.block(0, start, r, cols) = tmp2;
   }

   void IWorlandBackend::applyTriSolve(Matrix& out, const int start, const int cols, Matrix& A) const
   {
      char UPLO;
      if(A(0,0) == UPPER_BANDED)
      {
         UPLO = 'U';
      } else
      {
         assert(A(1,A.cols()-1) == LOWER_BANDED);
         UPLO = 'L';
      }
      char TRANS = 'N';
      char DIAG = 'N';
      int K = 1;
      int LDA = 2;
      int INCX = 1;
      int r = A.cols();
      for(int j = 0; j < cols; j++)
      {
         double *X = out.data() + (start+j)*out.rows();
         BLASFUNC(dtbsv)(&UPLO, &TRANS, &DIAG, &r, &K, A.data(), &LDA, X, &INCX);
      }
   }

   void IWorlandBackend::applyPair(Matrix& out, const int start, const int cols, const int rows, const SparseMatrix& P, const SparseMatrix& S) const
   {
      int mr = P.rows();
      int mc = P.cols();

      Matrix tmp = P*out.block(0, start, mc, cols);

      Eigen::SparseLU<SparseMatrix> solver;
      solver.compute(S);
      Matrix tmp2 = solver.solve(tmp);
      out.block(0, start, mr, cols) = tmp2;
   }

#if QUICC_IWORLANDBACKEND==QUICC_IWORLANDBACKEND_BLAS
   void IWorlandBackend::applyPair(Matrix& out, const int start, const int cols, const int rows, const Matrix& PS) const
   {
      int r = rows;
      char P_UPLO;
      char S_UPLO;
      if(PS(0,0) == UPPER_BANDED)
      {
         assert(PS(3,PS.cols()-1) == LOWER_BANDED);
         P_UPLO = 'U';
         S_UPLO = 'L';
      } else
      {
         assert(PS(1,PS.cols()-1) == LOWER_BANDED || PS(1,PS.cols()-1) == LOWER_BANDED_PADDED);
         assert(PS(2,0) == UPPER_BANDED);
         P_UPLO = 'L';
         S_UPLO = 'U';
         if(PS(1,PS.cols()-1) == LOWER_BANDED_PADDED)
         {
            out.block(r, start, 1, cols).setZero();
         }
      }
      char TRANS = 'N';
      char DIAG = 'N';
      int K = 1;
      int LDA = PS.rows();
      int INCX = 1;

      double *pP = const_cast<Matrix&>(PS).data();
      double *pS = pP+2;
      for(int j = 0; j < cols; j++)
      {
         double *X = out.data() + (start+j)*out.rows();
         BLASFUNC(dtbmv)(&P_UPLO, &TRANS, &DIAG, &r, &K, pP, &LDA, X, &INCX);
         BLASFUNC(dtbsv)(&S_UPLO, &TRANS, &DIAG, &r, &K, pS, &LDA, X, &INCX);
      }
   }

#else

   void IWorlandBackend::applyPair(Matrix& out, const int start, const int cols, const int rows, const Matrix& PS) const
   {
      if(PS(2,0) == UPPER_BANDED)
      {
#if QUICC_IWORLANDBACKEND==QUICC_IWORLANDBACKEND_SCALAR
         int dR = out.rows();
         MHDFloat *X = out.data() + start*dR;
         int q0 = 0;

         for(int j = 0;j < cols; ++j)
         {
            q0 = j*dR;
            applyUpperLower(X, q0, rows, PS);
         }
#elif QUICC_IWORLANDBACKEND==QUICC_IWORLANDBACKEND_VECTOR
         applyUpperLower(out.block(0, start, out.rows(), cols), rows, PS);
#endif

      } else
      {
         assert(PS(2,PS.cols()-1) == LOWER_BANDED || PS(2,PS.cols()-1) == LOWER_BANDED_PADDED);
         if(PS(2,PS.cols()-1) == LOWER_BANDED_PADDED)
         {
            out.block(rows, start, 1, cols).setZero();
         }
#if QUICC_IWORLANDBACKEND==QUICC_IWORLANDBACKEND_SCALAR
         int dR = out.rows();
         MHDFloat *X = out.data() + start*dR;
         int q0 = 0;

         for(int j = 0;j < cols; ++j)
         {
            q0 = j*dR + rows - 1;
            applyLowerUpper(X, q0, rows, PS);
         }
#elif QUICC_IWORLANDBACKEND==QUICC_IWORLANDBACKEND_VECTOR

         applyLowerUpper(out.block(0, start, out.rows(), cols), rows, PS);
#endif
      }
   }
#endif

   void IWorlandBackend::io(MHDFloat* out, const MHDFloat* in) const
   {
      this->mpOut = out;
      this->mpIn = in;
   }

   void IWorlandBackend::scaleC(const MHDFloat c, const bool isEven, const unsigned int id) const
   {
      Matrix& wTmp = this->workTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         wTmp.block(0, start, rows, cols) *= c;
         start += cols;
      }
   }

   void IWorlandBackend::scaleALPY(const MHDFloat a, const MHDFloat y, const bool isEven, const int lshift, const unsigned int id) const
   {
      Matrix& wTmp = this->workTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven, id))
      {
         int l = std::get<4>(loc) + lshift;
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         wTmp.block(0, start, rows, cols) *= (a*l + y);
         start += cols;
      }
   }

   void IWorlandBackend::scaleD(const bool isEven, const int lshift, const unsigned int id) const
   {
      const MHDFloat alpha = -0.5;
      Matrix& wTmp = this->workTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven, id))
      {
         int l = std::get<4>(loc) + lshift;
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         for(int n = 0; n < rows; n++)
         {
            wTmp.block(n, start, 1, cols) *= 2.0*(l + n + 1.0)*std::sqrt((n+1.0)/(n + alpha + l + 1.5));
         }
         start += cols;
      }
   }

   void IWorlandBackend::scaleSphLaplA(const bool isEven, const int lshift, const unsigned int id) const
   {
      const MHDFloat alpha = -0.5;
      Matrix& wTmp = this->workTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven, id))
      {
         int l = std::get<4>(loc) + lshift;
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         for(int n = 0; n < rows; n++)
         {
            wTmp.block(n, start, 1, cols) *= 4.0*(l + n + 2.0)*(l + n + 3.0)*std::sqrt((n+1.0)*(n+2.0)/((n + alpha + l + 2.5)*(n + alpha + l + 3.5)));
         }
         start += cols;
      }
   }

   void IWorlandBackend::scaleSphLaplB(const bool isEven, const int lshift, const unsigned int id) const
   {
      const MHDFloat alpha = -0.5;
      Matrix& wTmp = this->workTmp(id);
      int start = 0;
      for(auto loc: *this->pLoc(isEven, id))
      {
         int l = std::get<4>(loc) + lshift;
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc);
         for(int n = 0; n < rows; n++)
         {
            wTmp.block(n, start, 1, cols) *= 2.0*(l + n + 1.0)*(2.0*l + 3.0)*std::sqrt((n+1.0)/(n + alpha + l + 1.5));
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
      Matrix& wTmp = this->workTmp(id);
      int s = std::abs(nshift);

      int start = 0;
      for(auto loc: *this->pLoc(isEven,id))
      {
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc) - s;
         if(nshift > 0)
         {
            Matrix tmp = wTmp.block(0, start, rows, cols);
            wTmp.block(s, start, rows, cols) = tmp;
            wTmp.topRows(s).setZero();
         } else
         {
            Matrix tmp = wTmp.block(s, start, rows, cols);
            wTmp.block(0, start, rows, cols) = tmp;
            wTmp.block(rows, start, s, cols).setZero();
         }
         start += cols;
      }
   }

   void IWorlandBackend::copy(const int to, const int from, const int nshift, const bool isEven) const
   {
      assert(to != from);
      Matrix& wTmp = this->workTmp(from);
      Matrix& extraTmp = this->workTmp(to);
      int s = std::abs(nshift);

      if(to > 0)
      {
         *this->pLoc(isEven,to) = *this->pLoc(isEven,from);
      }

      int start = 0;
      for(auto loc: *this->pLoc(isEven,to))
      {
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc) - s;
         if(nshift > 0)
         {
            extraTmp.block(s, start, rows, cols) = wTmp.block(0, start, rows, cols);
            extraTmp.topRows(s).setZero();
         } else
         {
            extraTmp.block(0, start, rows, cols) = wTmp.block(s, start, rows, cols);
            extraTmp.block(rows, start, s, cols).setZero();
         }
         start += cols;
      }
   }

   void IWorlandBackend::add(const int to, const int from, const int nshift, const bool isEven) const
   {
      Matrix& wTmp = this->workTmp(to);
      Matrix& extraTmp = this->workTmp(from);
      int s = std::abs(nshift);

      int start = 0;
      for(auto loc: *this->pLoc(isEven,to))
      {
         int cols = std::get<2>(loc);
         int rows = std::get<3>(loc) - s;
         if(nshift > 0)
         {
            wTmp.block(s, start, rows, cols) += extraTmp.block(0, start, rows, cols);
         } else
         {
            wTmp.block(0, start, rows, cols) += extraTmp.block(s, start, rows, cols);
         }
         start += cols;
      }

   }

   void applyUpperLower(double* X, int q0, const int r, const Matrix& PS)
   {
      X[q0] = (PS(0,0)*X[q0] + PS(1,0)*X[q0+1]);
      q0++;

      for(int i = 1; i < r-1; ++i,++q0)
      {
         X[q0] = (PS(0,i)*X[q0] + PS(1,i)*X[q0+1] + PS(2, i)*X[q0-1]);
      }

      X[q0] = (X[q0]*PS(0,r-1) + PS(2, r-1)*X[q0-1]);
   }

   void applyLowerUpper(double* X, int q0, const int r, const Matrix& PS)
   {
      X[q0] = (PS(0,r-1)*X[q0-1] + PS(1,r-1)*X[q0]);
      --q0;

      for(int i = r-2; i > 0; --i,--q0)
      {
         X[q0] = (PS(0,i)*X[q0-1] + PS(1,i)*X[q0] + PS(2,i)*X[q0+1]);
      }

      X[q0] = (PS(1,0)*X[q0] + PS(2,0)*X[q0+1]);
   }

   void applyUpperLower(Eigen::Ref<Matrix> out, const int r, const Matrix& PS)
   {
      out.row(0) = (PS(0,0)*out.row(0) + PS(1,0)*out.row(1));

      for(int i = 1; i < r-1; ++i)
      {
         out.row(i) = (PS(0,i)*out.row(i) + PS(1,i)*out.row(i+1) + PS(2, i)*out.row(i-1));
      }
      out.row(r-1) = (PS(0,r-1)*out.row(r-1) + PS(2, r-1)*out.row(r-2));
   }

   void applyLowerUpper(Eigen::Ref<Matrix> out, const int r, const Matrix& PS)
   {
      out.row(r-1) = (PS(0,r-1)*out.row(r-2) + PS(1,r-1)*out.row(r-1));

      for(int i = r-2; i > 0; --i)
      {
         out.row(i) = (PS(0,i)*out.row(i-1) + PS(1,i)*out.row(i) + PS(2,i)*out.row(i+1));
      }

      out.row(0) = (PS(1,0)*out.row(0) + PS(2,0)*out.row(1));
   }

}
}
}
}
}
