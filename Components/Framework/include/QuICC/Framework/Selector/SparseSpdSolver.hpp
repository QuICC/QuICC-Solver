/** 
 * @file SparseSpdSolver.hpp
 * @brief Selector for configured sparse SPD solver
 */

#ifndef QUICC_TRANSFORM_SELECTOR_SPARSESPDSOLVER_HPP
#define QUICC_TRANSFORM_SELECTOR_SPARSESPDSOLVER_HPP

#if defined QUICC_SPSPDLINALG_UMFPACK
#include "QuICC/Solver/Sparse/UmfPackLU.hpp"
#elif defined QUICC_SPSPDLINALG_SPARSELU
#include "QuICC/Solver/Sparse/SparseLUCOLAMD.hpp"
#elif defined QUICC_SPSPDLINALG_MUMPS
#include "QuICC/Solver/Sparse/MumpsLU.hpp"
#elif defined QUICC_SPSPDLINALG_SIMPLICIALLLT
#include "QuICC/Solver/Sparse/SimplicialLLT.hpp"
#elif defined QUICC_SPSPDLINALG_SIMPLICIALLDLT
#include "QuICC/Solver/Sparse/SimplicialLDLT.hpp"
#endif //defined QUICC_SPSPDLINALG_UMFPACK

namespace QuICC {

namespace Transform {

namespace Selector {

   #if defined QUICC_SPSPDLINALG_UMFPACK
   template <class TMatrix>
      using SparseSpdSolver = QuICC::Solver::Sparse::UmfPackLU<TMatrix>;
   #elif defined QUICC_SPSPDLINALG_SPARSELU
   template <class TMatrix>
      using SparseSpdSolver = QuICC::Solver::Sparse::SparseLUCOLAMD<TMatrix>;
   #elif defined QUICC_SPSPDLINALG_MUMPS
   template <class TMatrix>
      using SparseSpdSolver = QuICC::Solver::Sparse::MumpsLU<TMatrix>;
   #elif defined QUICC_SPSPDLINALG_SIMPLICIALLT
   template <class TMatrix>
      using SparseSpdSolver = QuICC::Solver::Sparse::SimplicialLLT<TMatrix>;
   #elif defined QUICC_SPSPDLINALG_SIMPLICIALDLT
   template <class TMatrix>
      using SparseSpdSolver = QuICC::Solver::Sparse::SimplicialLDLT<TMatrix>;
   #endif //defined QUICC_SPSPDLINALG_UMFPACK

}
}
}

#endif // QUICC_TRANSFORM_SELECTOR_SPARSESPDSOLVER_HPP
