/** 
 * @file DifferentialSolver.hpp
 * @brief Interface for a differentiate by linear solver
 */

#ifndef QUICC_TRANSFORM_FFT_BACKEND_FFTW_DIFFERENTIALSOLVER_HPP
#define QUICC_TRANSFORM_FFT_BACKEND_FFTW_DIFFERENTIALSOLVER_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

   /**
    * @brief Differential by linear solver
    */
   class DifferentialSolver {
      public:
         /**
          * @brief Constructor
          */
         DifferentialSolver(const int specSize, const int blockSize, const int extraRows = 0);

         /**
          * @brief Destructor
          */
         ~DifferentialSolver();

         /**
          * @brief Set solver input
          */
         void input(const Matrix& in, const int shift = 0);

         /**
          * @brief Set solver input
          */
         void input(const MatrixZ& in, const int shift, const bool useReal);

         /**
          * @brief Set solver input multiplied by spectral operator
          */
         void inputSpectral(const Matrix& in);

         /**
          * @brief Set solver input multiplied by spectral operator
          */
         void inputSpectral(const int shift);

         /**
          * @brief Set solver input multiplied by spectral operator
          */
         void inputSpectral(const MatrixZ& in, const bool useReal);

         /**
          * @brief Set solver operator
          */
         void setOperator(const SparseMatrix& mat, const int extraRows = 0, const int extraCols = 0);

         /**
          * @brief Set solver operator
          */
         void setSpectralOperator(const SparseMatrix& mat, const int shift = 0);

         /**
          * @brief Set solver operator
          */
         void solve(const int zeroRows, Matrix& rOut);

      private:
         /**
          * @brief Spec size
          */
         int mSpecSize;

         /**
          * @brief Block size
          */
         int mBlockSize;

         /**
          * @brief Storage for solver input/output
          */
         Matrix mTmp;

         /**
          * @brief Sparse matrix used in solver
          */
         SparseMatrix mOp;

         /**
          * @brief Spectral operator
          */
         SparseMatrix mSpecOp;

         /// Forward declaration for solver
         class SolverImpl;

         /**
          * @brief Sparse linear solver
          */
         std::unique_ptr<SolverImpl> mpSolverImpl;
   };

}
}
}
}
}

#endif // QUICC_TRANSFORM_FFT_BACKEND_FFTW_DIFFERENTIALSOLVER_HPP
