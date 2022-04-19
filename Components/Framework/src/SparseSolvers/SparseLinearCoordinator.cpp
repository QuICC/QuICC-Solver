/**
 * @file SparseLinearCoordinator.cpp
 * @brief Implementation of the sparse linear coordinator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SparseSolvers/SparseLinearCoordinator.hpp"

// Project includes
//

namespace QuICC {

namespace Solver {

   SparseLinearCoordinator::SparseLinearCoordinator()
      : SparseLinearCoordinatorBase<SparseLinearSolver>()
   {
   }

   SparseLinearCoordinator::~SparseLinearCoordinator()
   {
   }

   void SparseLinearCoordinator::solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq, scalVar, vectVar);

      // Solve all the linear systems
      this->solveSystems();

      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);

      // Clear the solver RHS
      this->clearSolvers();
   }

   void SparseLinearCoordinator::buildSolverMatrix(SharedRealSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildLinearSolverMatrixWrapper(spSolver, spEq, comp, idx);
   }

   void SparseLinearCoordinator::buildSolverMatrix(SharedComplexSolverType spSolver, Equations::SharedIEquation spEq, FieldComponents::Spectral::Id comp, const int idx)
   {
      buildLinearSolverMatrixWrapper(spSolver, spEq, comp, idx);
   }

}
}
