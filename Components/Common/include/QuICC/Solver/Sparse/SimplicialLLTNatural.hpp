/** 
 * @file SimplicialLLTNatural.hpp
 * @brief Simple typedef to use SimplicialLLTNatural solver.
 */

#ifndef QUICC_SOLVER_SPARSE_SIMPLICIALLLTNATURAL_HPP
#define QUICC_SOLVER_SPARSE_SIMPLICIALLLTNATURAL_HPP

#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using SimplicialLLTNatural = Eigen::SimplicialLLT<TMatrix, Eigen::Lower, Eigen::NaturalOrdering<int> >;
}
}
}

#endif // QUICC_SOLVER_SPARSE_SIMPLICIALLLTNATURAL_HPP
