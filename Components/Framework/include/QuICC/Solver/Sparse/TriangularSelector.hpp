/** 
 * @file TriangularSelector.hpp
 * @brief Small template to select different sparse triangular solvers.
 */

#ifndef QUICC_SOLVER_SPARSE_TRIANGULARSELECTOR_HPP
#define QUICC_SOLVER_SPARSE_TRIANGULARSELECTOR_HPP

// UmfPack Version for triangular solve
#ifdef QUICC_SPTRILINALG_UMFPACK

#include <Eigen/UmfPackSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use UmfPack as Triangular sparse solver 
    */
   template <typename TMatrix> struct TriangularSelector
   {
      typedef Eigen::UmfPackLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPTRILINALG_UMFPACK

// SparseLU Version for triangular solve
#ifdef QUICC_SPTRILINALG_SPARSELU

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use SparseLU as Triangular sparse solver 
    */
   template <typename TMatrix> struct TriangularSelector
   {
      typedef Eigen::SparseLU<TMatrix, Eigen::COLAMDOrdering<int> > Type;
   };

}
}
}
#endif //QUICC_SPTRILINALG_SPARSELU

// MUMPS Version for triangular solve
#ifdef QUICC_SPTRILINALG_MUMPS

#include "../External/Interfaces/MumpsLU.hpp"

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use MUMPS as Triangular sparse solver 
    */
   template<typename TMatrix> struct TriangularSelector
   {
      typedef Eigen::MumpsLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPTRILINALG_MUMPS

#endif // QUICC_SOLVER_SPARSE_TRIANGULARSELECTOR_HPP
