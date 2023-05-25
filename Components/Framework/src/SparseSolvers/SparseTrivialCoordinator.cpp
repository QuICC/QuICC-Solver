/**
 * @file SparseTrivialCoordinator.cpp
 * @brief Implementation of the sparse trivial coordinator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSolvers/SparseTrivialCoordinator.hpp"

namespace QuICC {

namespace Solver {

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

      //
      // Initialise the solver storage
      //

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

      //
      // Initialise the solvers initial state
      //

      // Initialise with initial state
      this->initSolution(scalEq, vectEq);
   }

} // Solver
} // QuICC
