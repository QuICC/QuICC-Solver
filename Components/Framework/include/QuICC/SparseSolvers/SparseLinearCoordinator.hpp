/**
 * @file SparseLinearCoordinator.hpp
 * @brief Implementation of a general sparse linear solver coordinator
 */

#ifndef QUICC_SOLVER_SPARSELINEARCOORDINATOR_HPP
#define QUICC_SOLVER_SPARSELINEARCOORDINATOR_HPP

// Configuration includes
//

// Debug includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/SparseSolvers/SparseLinearCoordinatorBase.hpp"
#include "QuICC/SparseSolvers/SparseLinearSolver.hpp"
#include "QuICC/Equations/IScalarEquation.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Solver {

   /**
    * @brief Implementation of general sparse linear solver coordinator
    */
   class SparseLinearCoordinator: public SparseLinearCoordinatorBase<SparseLinearSolver>
   {
      public:
         /// Typedef for a shared scalar equation range
         typedef typename SparseLinearCoordinatorBase<SparseLinearSolver>::ScalarEquation_range  ScalarEquation_range;

         /// Typedef for a shared vector equation range
         typedef typename SparseLinearCoordinatorBase<SparseLinearSolver>::VectorEquation_range  VectorEquation_range;

         /// Typedef for a shared scalar variable map
         typedef typename SparseLinearCoordinatorBase<SparseLinearSolver>::ScalarVariable_map  ScalarVariable_map;

         /// Typedef for a shared vector variable map
         typedef typename SparseLinearCoordinatorBase<SparseLinearSolver>::VectorVariable_map  VectorVariable_map;

         /// Typedef for a shared real operator solver
         typedef typename SparseLinearCoordinatorBase<SparseLinearSolver>::SharedRealSolverType  SharedRealSolverType;

         /// Typedef for a shared complex operator solver
         typedef typename SparseLinearCoordinatorBase<SparseLinearSolver>::SharedComplexSolverType  SharedComplexSolverType;

         /**
          * @brief Constructor
          */
         SparseLinearCoordinator();

         /**
          * @brief Destructor
          */
         ~SparseLinearCoordinator();

         /**
          * @brief Solve the equations
          *
          * @param scalEq Shared scalar equations
          * @param vectEq Shared vector equations
          * @param scalVar Shared scalar variables
          * @param vectVar Shared vector variables
          */
         void solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar);

      protected:
         /**
          * @brief Build the real operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

         /**
          * @brief Build the complex operator, complex field solver matrix
          *
          * @param spSolver   Shared sparse real solver
          * @param spEq       Shared pointer to equation
          * @param comp       Field component
          * @param idx        Matrix index
          */
         virtual void buildSolverMatrix(SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

      private:
   };

   /**
    * @brief Generic implementation to build solver matrices
    */
   template <typename TSolver> void buildLinearSolverMatrixWrapper(std::shared_ptr<TSolver > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx);

   //
   //
   //

   template <typename TSolver> void buildLinearSolverMatrixWrapper(std::shared_ptr<TSolver > spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      std::map<std::size_t, DecoupledZSparse> ops;

      // Get model's linear operator with tau lines
      ops.insert(std::make_pair(ModelOperator::ImplicitLinear::id(), DecoupledZSparse()));
      spEq->buildModelMatrix(ops.find(ModelOperator::ImplicitLinear::id())->second, ModelOperator::ImplicitLinear::id(), comp, idx, ModelOperatorBoundary::SolverHasBc::id());

      spSolver->buildOperators(idx, ops, spEq->couplingInfo(comp).systemN(idx));

      // Solver is initialized
      spSolver->setInitialized();
   }
}
}

#endif // QUICC_SOLVER_SPARSELINEARCOORDINATOR_HPP
