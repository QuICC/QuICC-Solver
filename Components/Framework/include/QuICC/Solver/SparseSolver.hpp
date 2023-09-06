/** 
 * @file SparseSolver.hpp
 * @brief Selector for configured sparse solver
 */

#ifndef QUICC_FRAMEWORK_SELECTOR_SPARSESOLVER_HPP
#define QUICC_FRAMEWORK_SELECTOR_SPARSESOLVER_HPP

#if defined QUICC_SPLINALG_SUPERLU
#include "QuICC/Solver/Sparse/SuperLU.hpp"
#elif defined QUICC_SPLINALG_UMFPACK
#include "QuICC/Solver/Sparse/UmfPackLU.hpp"
#elif defined QUICC_SPLINALG_SPARSELUNATURAL
#include "QuICC/Solver/Sparse/SparseLUNATURAL.hpp"
#elif defined QUICC_SPLINALG_SPARSELUAMD
#include "QuICC/Solver/Sparse/SparseLUAMD.hpp"
#elif defined QUICC_SPLINALG_SPARSELU || defined QUICC_SPLINALG_SPARSELUCOLAMD
#include "QuICC/Solver/Sparse/SparseLUCOLAMD.hpp"
#elif defined QUICC_SPLINALG_SPARSELUMETIS
#include "QuICC/Solver/Sparse/SparseLUMetis.hpp"
#elif defined QUICC_SPLINALG_MUMPS
#include "QuICC/Solver/Sparse/MumpsLU.hpp"
#endif //defined QUICC_SPLINALG_SPARSELU

namespace QuICC {

namespace Framework {

namespace Selector {

   #if defined QUICC_SPLINALG_SUPERLU
   template <class TMatrix>
      using SparseSolver = QuICC::Solver::Sparse::SuperLU<TMatrix>;
   #elif defined QUICC_SPLINALG_UMFPACK
   template <class TMatrix>
      using SparseSolver = QuICC::Solver::Sparse::UmfPackLU<TMatrix>;
   #elif defined QUICC_SPLINALG_SPARSELUNATURAL
   template <class TMatrix>
      using SparseSolver = QuICC::Solver::Sparse::SparseLUNATURAL<TMatrix>;
   #elif defined QUICC_SPLINALG_SPARSELUAMD
   template <class TMatrix>
      using SparseSolver = QuICC::Solver::Sparse::SparseLUAMD<TMatrix>;
   #elif defined QUICC_SPLINALG_SPARSELUCOLAMD || defined QUICC_SPLINALG_SPARSELU
   template <class TMatrix>
      using SparseSolver = QuICC::Solver::Sparse::SparseLUCOLAMD<TMatrix>;
   #elif defined QUICC_SPLINALG_SPARSELUMETIS
   template <class TMatrix>
      using SparseSolver = QuICC::Solver::Sparse::SparseLUMetis<TMatrix>;
   #elif defined QUICC_SPLINALG_MUMPS
   template <class TMatrix>
      using SparseSolver = QuICC::Solver::Sparse::MumpsLU<TMatrix>;
   #endif //defined QUICC_SPLINALG_SPARSELU

}
}
}

#endif // QUICC_FRAMEWORK_SELECTOR_SPARSESOLVER_HPP
