/** 
 * @file SparseTriangularSolver.hpp
 * @brief Selector for configured triangular solver used with FFTW backend
 */

#ifndef QUICC_TRANSFORM_SELECTOR_SPARSETRIANGULARSOLVER_HPP
#define QUICC_TRANSFORM_SELECTOR_SPARSETRIANGULARSOLVER_HPP

#if defined QUICC_TRANSFORM_FFTW_SPTRILINALG_UMFPACK
#include "QuICC/Solver/Sparse/UmfPackLU.hpp"
#elif defined QUICC_TRANSFORM_FFTW_SPTRILINALG_SPARSELU
#include "QuICC/Solver/Sparse/SparseLUCOLAMD.hpp"
#elif defined QUICC_TRANSFORM_FFTW_SPTRILINALG_MUMPS
#include "QuICC/Solver/Sparse/MumpsLU.hpp"
#endif //defined QUICC_TRANSFORM_FFTW_SPTRILINALG_UMFPACK

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Backend {

namespace Fftw {

namespace Selector {

   #if defined QUICC_TRANSFORM_FFTW_SPTRILINALG_UMFPACK
   template <class TMatrix>
      using SparseTriangularSolver = QuICC::Solver::Sparse:UmfPackLU<TMatrix>;
   #elif defined QUICC_TRANSFORM_FFTW_SPTRILINALG_SPARSELU
   template <class TMatrix>
      using SparseTriangularSolver = QuICC::Solver::Sparse::SparseLUCOLAMD<TMatrix>;
   #elif defined QUICC_TRANSFORM_FFTW_SPTRILINALG_MUMPS
   template <class TMatrix>
      using SparseTriangularSolver = QuICC::Solver::Sparse:MumpsLU<TMatrix>;
   #endif //defined QUICC_TRANSFORM_FFTW_SPTRILINALG_UMFPACK

}
}
}
}
}
}

#endif // QUICC_TRANSFORM_SELECTOR_SPARSETRIANGULARSOLVER_HPP
