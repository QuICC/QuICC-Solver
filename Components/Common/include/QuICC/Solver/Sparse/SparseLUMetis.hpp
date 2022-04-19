/** 
 * @file SparseLUMetis.hpp
 * @brief Simple typedef to use SparseLUMetis solver.
 */

#ifndef QUICC_SOLVER_SPARSE_SPARSELUMETIS_HPP
#define QUICC_SOLVER_SPARSE_SPARSELUMETIS_HPP

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>
#include <Eigen/MetisSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using SparseLUMetis = Eigen::SparseLU<TMatrix, Eigen::MetisOrdering<int> >;
}
}
}

#endif // QUICC_SOLVER_SPARSE_SPARSELUMETIS_HPP
