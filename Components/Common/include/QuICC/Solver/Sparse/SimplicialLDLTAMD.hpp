/** 
 * @file SimplicialLLTAMD.hpp
 * @brief Simple typedef to use SimplicialLLTAMD solver.
 */

#ifndef QUICC_SOLVER_SPARSE_SIMPLICIALLLTAMD_HPP
#define QUICC_SOLVER_SPARSE_SIMPLICIALLLTAMD_HPP

#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using SimplicialLLTAMD = Eigen::SimplicialLLT<TMatrix, Eigen::Lower, Eigen::AMDOrdering<int> >;
}
}
}

#endif // QUICC_SOLVER_SPARSE_SIMPLICIALLLTAMD_HPP
