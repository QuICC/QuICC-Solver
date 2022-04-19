/** 
 * @file GeneralSelector.hpp
 * @brief Small template to select different sparse general solvers.
 */

#ifndef QUICC_SOLVER_SPARSE_GENERALSELECTOR_HPP
#define QUICC_SOLVER_SPARSE_GENERALSELECTOR_HPP

// SuperLU Version for general solve
#ifdef QUICC_SPLINALG_SUPERLU

#include <Eigen/SuperLUSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use SuperLU as sparse general solver 
    */
   template <typename TMatrix> struct GeneralSelector
   {
      typedef Eigen::SuperLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPLINALG_SUPERLU

// UmfPack Version for general solve
#ifdef QUICC_SPLINALG_UMFPACK

#include <Eigen/UmfPackSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use UmfPack as sparse general solver 
    */
   template <typename TMatrix> struct GeneralSelector
   {
      typedef Eigen::UmfPackLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPLINALG_UMFPACK

// SparseLU Version for general solve
#ifdef QUICC_SPLINALG_SPARSELU

#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>
//#include <Eigen/MetisSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use Eigen's SparseLU as sparse general solver
    */
   template <typename TMatrix> struct GeneralSelector
   {
      //typedef Eigen::SparseLU<TMatrix, Eigen::NaturalOrdering<int> > Type;
      //typedef Eigen::SparseLU<TMatrix, Eigen::AMDOrdering<int> > Type;
      typedef Eigen::SparseLU<TMatrix, Eigen::COLAMDOrdering<int> > Type;
      //typedef Eigen::SparseLU<TMatrix, Eigen::MetisOrdering<int> > Type;
   };

}
}
}
#endif //QUICC_SPLINALG_SPARSELU

// KLU Version for general solve
#ifdef QUICC_SPLINALG_KENTLU

#warning "THIS SOLVER NEED EVALUTION WETHER IT IS USEFUL"
#include <Eigen/KentLUSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use KentLU as sparse general solver 
    */
   template <typename TMatrix> struct GeneralSelector
   {
      typedef Eigen::KentLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPLINALG_KENTLU

// SparseQR Version for general solve
#ifdef QUICC_SPLINALG_SPARSEQR

#warning "THIS SOLVER NEED EVALUTION WETHER IT IS USEFUL"
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
//#include <Eigen/MetisSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use SparseQR as sparse general solver
    */
   template<typename TMatrix> struct GeneralSelector
   {
      //typedef Eigen::SparseQR<TMatrix, Eigen::NaturalOrdering<int> > Type;
      //typedef Eigen::SparseQR<TMatrix, Eigen::AMDOrdering<int> > Type;
      typedef Eigen::SparseQR<TMatrix, Eigen::COLAMDOrdering<int> > Type;
      //typedef Eigen::SparseQR<TMatrix, Eigen::MetisOrdering<int> > Type;
   };

}
}
}
#endif //QUICC_SPLINALG_SPARSEQR

// SuiteSparseQR Version for general solve
#ifdef QUICC_SPLINALG_SPQR

#warning "THIS SOLVER NEED EVALUTION WETHER IT IS USEFUL"
#include <Eigen/SPQRSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use SuiteSparseQR (SPQR) as sparse general solver 
    */
   template<typename TMatrix> struct GeneralSelector
   {
      typedef Eigen::SPQR<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPLINALG_SPQR

// Pardiso Version for general solve
#ifdef QUICC_SPLINALG_PARDISO

#warning "THIS SOLVER NEED EVALUTION WETHER IT IS USEFUL"
#include "../External/Interfaces/PardisoLU.hpp"

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use Pardiso as general sparse Solver 
    */
   template<typename TMatrix> struct GeneralSelector
   {
      typedef Eigen::PardisoLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPLINALG_PARDISO

// MKL Pardiso Version for general solve
#ifdef QUICC_SPLINALG_MKLPARDISO

#warning "THIS SOLVER NEED EVALUTION WETHER IT IS USEFUL"
#include <Eigen/PardisoSupport>

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use MKL's Pardiso as general sparse solver 
    */
   template<typename TMatrix> struct GeneralSelector
   {
      typedef Eigen::PardisoLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPLINALG_MKLPARDISO

// MUMPS Version for general solve
#ifdef QUICC_SPLINALG_MUMPS

#include "../External/Interfaces/MumpsLU.hpp"

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use MUMPS as sparse general solver 
    */
   template<typename TMatrix> struct GeneralSelector
   {
      typedef Eigen::MumpsLU<TMatrix> Type;
   };

}
}
}
#endif //QUICC_SPLINALG_MUMPS

// BiCGSTAB Version for general solve
#ifdef QUICC_SPLINALG_BICGSTAB

#warning "THIS SOLVER NEED EVALUTION WETHER IT IS USEFUL"
#include <Eigen/IterativeLinearSolvers>   

namespace QuICC {

namespace Solver {

namespace Sparse {

   /**
    * @brief Use Eigen's BiCGSTAB as general sparse solver
    */
   template<typename TMatrix> struct GeneralSelector
   {
      //typedef Eigen::BiCGSTAB<TMatrix, Eigen::IdentityPreconditioner<typename TMatrix::Scalar> > Type;
      typedef Eigen::BiCGSTAB<TMatrix, Eigen::DiagonalPreconditioner<typename TMatrix::Scalar> > Type;
      //typedef Eigen::BiCGSTAB<TMatrix, Eigen::IncompleteLUT<typename TMatrix::Scalar> > Type;
   };

}
}
}
#endif //QUICC_SPLINALG_BICGSTAB

#endif // QUICC_SOLVER_SPARSE_GENERALSELECTOR_HPP
