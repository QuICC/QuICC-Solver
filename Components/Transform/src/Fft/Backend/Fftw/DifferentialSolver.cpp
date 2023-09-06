/**
 * @file ChebyshevProjector.cpp
 * @brief Source of the interface to differential by linear solver
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Transform/Fft/Backend/Fftw/DifferentialSolver.hpp"

// Project includes
//
#include "QuICC/Transform/Fft/Backend/Fftw/Selector/SparseTriangularSolver.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   class DifferentialSolver::SolverImpl {
      public:
         /// Typedef for solver
         typedef Selector::SparseTriangularSolver<SparseMatrix> SolverType;

         // Get solver
         SolverType& get() {return this->mSolver;};
      private:
         /// Solver
         SolverType mSolver;
   };

   DifferentialSolver::DifferentialSolver(const int specSize, const int blockSize, const int extraRows)
      : mSpecSize(specSize), mBlockSize(blockSize), mExtraRows(extraRows)
   {
      // Initializer solver storage
      this->mTmp.setZero(this->mSpecSize + extraRows, this->mBlockSize);

      this->mpSolverImpl = std::make_unique<SolverImpl>();
   }

   DifferentialSolver::~DifferentialSolver()
   {
   }

   void DifferentialSolver::input(const Matrix& in, const int shift)
   {
      this->mTmp.topRows(this->mSpecSize - shift) = in.block(shift, 0, this->mSpecSize - shift, in.cols());
   }

   void DifferentialSolver::inputSpectral(const Matrix& in)
   {
      this->mTmp.topRows(this->mSpecSize) = this->mSpecOp*in.topRows(this->mSpecSize);
   }

   void DifferentialSolver::inputSpectral(const int shift)
   {
      this->mTmp.topRows(this->mTmp.rows() - shift) = this->mSpecOp*this->mTmp;
   }

   void DifferentialSolver::input(const MatrixZ& in, const int shift, const bool useReal)
   {
      if(useReal)
      {
         this->mTmp.topRows(this->mSpecSize - shift) = in.real().block(shift, 0, this->mSpecSize - shift, in.cols());
      } else
      {
         this->mTmp.topRows(this->mSpecSize - shift) = in.imag().block(shift, 0, this->mSpecSize - shift, in.cols());
      }
   }

   void DifferentialSolver::inputSpectral(const MatrixZ& in, const bool useReal)
   {
      if(useReal)
      {
         this->mTmp.topRows(this->mSpecOp.rows()) = this->mSpecOp*in.real().topRows(this->mSpecSize);
      } else
      {
         this->mTmp.topRows(this->mSpecOp.rows()) = this->mSpecOp*in.imag().topRows(this->mSpecSize);
      }
   }

   void DifferentialSolver::setOperator(const SparseMatrix& mat, const int extraRows, const int extraCols)
   {
      this->mOp = mat.bottomLeftCorner(this->mSpecSize + extraRows, this->mSpecSize + extraCols);

      this->mpSolverImpl->get().compute(this->mOp);
      if(this->mpSolverImpl->get().info() != Eigen::Success)
      {
         throw std::logic_error("Factorization of transform operator failed!");
      }
   }

   void DifferentialSolver::setSpectralOperator(const SparseMatrix& mat, const int shift)
   {
      this->mSpecOp = mat.bottomRows(mat.rows() - shift);
   }

   SparseMatrix& DifferentialSolver::getSpectralOperator()
   {
      return this->mSpecOp;
   }

   void DifferentialSolver::solve(Matrix& tmp, const int zeroRows)
   {
      tmp.middleRows(this->mSpecSize + this->mExtraRows - zeroRows, zeroRows).setZero();
      Matrix x = this->mpSolverImpl->get().solve(tmp.topRows(this->mSpecSize + this->mExtraRows));
      tmp.topRows(this->mSpecSize + this->mExtraRows) = x;
   }

}
}
}
}
}
