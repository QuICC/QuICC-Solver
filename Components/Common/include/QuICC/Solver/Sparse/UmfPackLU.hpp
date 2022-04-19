/** 
 * @file UmfPackLU.hpp
 * @brief Simple typedef to use UmfPackLU solver.
 */

#ifndef QUICC_SOLVER_SPARSE_UMFPACKLU_HPP
#define QUICC_SOLVER_SPARSE_UMFPACKLU_HPP

#include <Eigen/UmfPackSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using UmfPackLU = Eigen::UmfPackLU<TMatrix>;
}
}
}

#endif // QUICC_SOLVER_SPARSE_UMFPACKLU_HPP
