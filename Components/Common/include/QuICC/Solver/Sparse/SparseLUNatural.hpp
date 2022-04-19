/** 
 * @file SparseLUNatural.hpp
 * @brief Simple typedef to use SparseLUNatural solver.
 */

#ifndef QUICC_SOLVER_SPARSE_SPARSELUNATURAL_HPP
#define QUICC_SOLVER_SPARSE_SPARSELUNATURAL_HPP

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using SparseLUNatural = Eigen::SparseLU<TMatrix, Eigen::NaturalOrdering<int> >;
}
}
}

#endif // QUICC_SOLVER_SPARSE_SPARSELUNATURAL_HPP
