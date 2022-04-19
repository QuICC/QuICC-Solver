/** 
 * @file SparseLUAMD.hpp
 * @brief Simple typedef to use SparseLUAMD solver.
 */

#ifndef QUICC_SOLVER_SPARSE_SPARSELUAMD_HPP
#define QUICC_SOLVER_SPARSE_SPARSELUAMD_HPP

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using SparseLUAMD = Eigen::SparseLU<TMatrix, Eigen::AMDOrdering<int> >;
}
}
}

#endif // QUICC_SOLVER_SPARSE_SPARSELUAMD_HPP
