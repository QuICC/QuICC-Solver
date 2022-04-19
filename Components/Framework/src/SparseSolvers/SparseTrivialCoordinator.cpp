/**
 * @file SparseTrivialCoordinator.cpp
 * @brief Implementation of the sparse trivial coordinator
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SparseSolvers/SparseTrivialCoordinator.hpp"

// Project includes
//

namespace QuICC {

namespace Solver {

   SparseTrivialCoordinator::SparseTrivialCoordinator()
      : SparseCoordinatorBase<SparseTrivialSolver>()
   {
   }

   SparseTrivialCoordinator::~SparseTrivialCoordinator()
   {
   }

   void SparseTrivialCoordinator::solve(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq, const ScalarVariable_map& scalVar, const VectorVariable_map& vectVar)
   {
      // Update the equation input to the timestepper
      this->getInput(scalEq, vectEq, scalVar, vectVar);

      // Transfer timestep output back to equations
      this->transferOutput(scalEq, vectEq);

      // Clear the solver RHS
      this->clearSolvers();
   }

   void SparseTrivialCoordinator::init(const ScalarEquation_range& scalEq, const VectorEquation_range& vectEq)
   {
      //
      // Create real/complex solvers
      //

      DebuggerMacro_start("Trivial: create solvers", 2);
      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         // Get type information for the solvers
         this->createSolver(scalEqIt, FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         // Get type information for the solvers
         for(auto& compIt: make_range(vectEqIt->spectralRange()))
         {
            this->createSolver(vectEqIt, compIt);
         }
      }
      DebuggerMacro_stop("Trivial: create solvers t = ", 2);

      //
      // Initialise the solver storage
      //

      DebuggerMacro_start("Trivial: create storage", 2);
      // Loop over all scalar equations
      for(auto& scalEqIt: make_range(scalEq))
      {
         // Create storage
         this->createStorage(scalEqIt, FieldComponents::Spectral::SCALAR);
      }

      // Loop over all vector equations
      for(auto& vectEqIt: make_range(vectEq))
      {
         // Create storage
         for(auto& compIt: make_range(vectEqIt->spectralRange()))
         {
            this->createStorage(vectEqIt, compIt);
         }
      }

      // Initialise the start rows
      this->initStartRow();
      DebuggerMacro_stop("Trivial: create storage t = ", 2);

      //
      // Initialise the solvers initial state
      //

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

}
}
