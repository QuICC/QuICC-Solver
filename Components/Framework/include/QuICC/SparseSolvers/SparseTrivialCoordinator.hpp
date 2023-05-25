/**
 * @file SparseTrivialCoordinator.hpp
 * @brief Implementation of the base for a general sparse trivial solver coordinator
 */

#ifndef QUICC_SOLVER_SPARSETRIVIALCOORDINATOR_HPP
#define QUICC_SOLVER_SPARSETRIVIALCOORDINATOR_HPP

// System includes
//

// Project includes
//
#include "QuICC/SparseSolvers/SparseCoordinatorBase.hpp"
#include "QuICC/SparseSolvers/SparseTrivialSolver.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of the base for a general sparse trivial solver coordinator
    */
   class SparseTrivialCoordinator: public SparseCoordinatorBase<SparseTrivialSolver>
   {
      public:
         /// Typedef for a shared scalar equation range
         typedef typename SparseCoordinatorBase<SparseTrivialSolver>::ScalarEquation_range  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef typename SparseCoordinatorBase<SparseTrivialSolver>::VectorEquation_range  VectorEquation_range;

         /// Typedef for a shared scalar variable map
         typedef typename SparseCoordinatorBase<SparseTrivialSolver>::ScalarVariable_map  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef typename SparseCoordinatorBase<SparseTrivialSolver>::VectorVariable_map  VectorVariable_map;

         /**
          * @brief Constructor
          */
         SparseTrivialCoordinator() = default;

         /**
          * @brief Destructor
          */
         virtual ~SparseTrivialCoordinator() = default;

         /**
          * @brief Initialise solver coordinator
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          */
         virtual void init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq);

         /**
          * @brief Solve the equations
          *
          * @param scalEq  Shared scalar equations
          * @param vectEq  Shared vector equations
          * @param scalVar Shared scalar variables
          * @param vectVar Shared vector variables
          */
         void solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

      protected:

      private:
   };
} // Solver
} // QuICC

#endif // QUICC_SOLVER_SPARSETRIVIALCOORDINATORBASE_HPP
