/** 
 * @file SparseLUCOLAMD.hpp
 * @brief Simple typedef to use SparseLUCOLAMD solver.
 */

#ifndef QUICC_SOLVER_SPARSE_SPARSELUCOLAMD_HPP
#define QUICC_SOLVER_SPARSE_SPARSELUCOLAMD_HPP

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

      template <class TMatrix>
         using SparseLUCOLAMD = Eigen::SparseLU<TMatrix, Eigen::COLAMDOrdering<int> >;
}
}
}

#endif // QUICC_SOLVER_SPARSE_SPARSELUCOLAMD_HPP
