/** 
 * @file SpdSelector.hpp
 * @brief Small template to select different sparse SPD solvers.
 */

#ifndef QUICC_SOLVER_SPARSE_SPDSELECTOR_HPP
#define QUICC_SOLVER_SPARSE_SPDSELECTOR_HPP

// UmfPack Version for SPD solve
#ifdef QUICC_SPSPDLINALG_UMFPACK

#include <Eigen/UmfPackSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use UmfPack as sparse SPD solver
    */
   template <typename TMatrix> struct SpdSelector
   {
      typedef Eigen::UmfPackLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPSPDLINALG_UMFPACK

// SparseLU Version for SPD solve
#ifdef QUICC_SPSPDLINALG_SPARSELU

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use Eigen's SparseLU as sparse SPD solver 
    */
   template <typename TMatrix> struct SpdSelector
   {
      typedef Eigen::SparseLU<TMatrix, Eigen::COLAMDOrdering<int> > Type;
   };

}
}
}
#endif //QUICC_SPSPDLINALG_SPARSELU

// MUMPS Version for SPD solve
#ifdef QUICC_SPSPDLINALG_MUMPS

#include "../External/Interfaces/MumpsLU.hpp"

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use MUMPS as sparse SPD solver 
    */
   template<typename TMatrix> struct SpdSelector
   {
      typedef Eigen::MumpsLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPSPDLINALG_MUMPS

// SparseLU Version for SPD solve
#ifdef QUICC_SPSPDLINALG_SIMPLICIALLLT

#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use Eigen's SimplicialLLT as sparse SPD solver 
    */
   template <typename TMatrix> struct SpdSelector
   {
      typedef Eigen::SimplicialLLT<TMatrix, Eigen::Lower, Eigen::AMDOrdering<int> > Type;
   };

}
}
}
#endif //QUICC_SPSPDLINALG_SIMPLICIALLLT

// SimplicialLDLT Version for SPD solve
#ifdef QUICC_SPSPDLINALG_SIMPLICIALLDLT

#include <Eigen/SparseCholesky>
#include <Eigen/OrderingMethods>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use Eigen's SimplicialLDLT as sparse SPD solver 
    */
   template <typename TMatrix> struct SpdSelector
   {
      typedef Eigen::SimplicialLDLT<TMatrix, Eigen::Lower, Eigen::AMDOrdering<int> > Type;
   };

}
}
}
#endif //QUICC_SPSPDLINALG_SIMPLICIALLDLT

#endif // QUICC_SOLVER_SPARSE_SPDSELECTOR_HPP
