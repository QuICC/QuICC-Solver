/** 
 * @file SuperLU.hpp
 * @brief Simple typedef to use SuperLU solver.
 */

#ifndef QUICC_SOLVER_SPARSE_SUPERLU_HPP
#define QUICC_SOLVER_SPARSE_SUPERLU_HPP

#include <Eigen/SuperLUSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using SuperLU = Eigen::SuperLU<TMatrix>;
}
}
}

#endif // QUICC_SOLVER_SPARSE_SUPERLU_HPP
