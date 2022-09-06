/** 
 * @file SparseTriangularSolver.hpp
 * @brief Selector for configured sparse triangular solver
 */

#ifndef QUICC_FRAMEWORK_SELECTOR_SPARSETRIANGULARSOLVER_HPP
#define QUICC_FRAMEWORK_SELECTOR_SPARSETRIANGULARSOLVER_HPP

#if defined QUICC_SPTRILINALG_UMFPACK
#include "QuICC/Solver/Sparse/UmfPackLU.hpp"
#elif defined QUICC_SPTRILINALG_SPARSELU
#include "QuICC/Solver/Sparse/SparseLUCOLAMD.hpp"
#elif defined QUICC_SPTRILINALG_MUMPS
#include "QuICC/Solver/Sparse/MumpsLU.hpp"
#endif //defined QUICC_SPTRILINALG_UMFPACK

namespace QuICC {

namespace Framework {

namespace Selector {

   #if defined QUICC_SPTRILINALG_UMFPACK
   template <class TMatrix>
      using SparseTriangularSolver = QuICC::Solver::Sparse::UmfPackLU<TMatrix>;
   #elif defined QUICC_SPTRILINALG_SPARSELU
   template <class TMatrix>
      using SparseTriangularSolver = QuICC::Solver::Sparse::SparseLUCOLAMD<TMatrix>;
   #elif defined QUICC_SPTRILINALG_MUMPS
   template <class TMatrix>
      using SparseTriangularSolver = QuICC::Solver::Sparse::MumpsLU<TMatrix>;
   #endif //defined QUICC_SPTRILINALG_UMFPACK

}
}
}

#endif // QUICC_FRAMEWORK_SELECTOR_SPARSETRIANGULARSOLVER_HPP
