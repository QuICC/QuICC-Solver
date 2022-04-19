/** 
 * @file SimplicialLDLTNatural.hpp
 * @brief Simple typedef to use SimplicialLDLTNatural solver.
 */

#ifndef QUICC_SOLVER_SPARSE_SIMPLICIALLDLTNATURAL_HPP
#define QUICC_SOLVER_SPARSE_SIMPLICIALLDLTNATURAL_HPP

#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using SimplicialLDLTNatural = Eigen::SimplicialLDLT<TMatrix, Eigen::Lower, Eigen::NaturalOrdering<int> >;
}
}
}

#endif // QUICC_SOLVER_SPARSE_SIMPLICIALLDLTNATURAL_HPP
