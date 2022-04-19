/**
 * @file StateGenerator.cpp
 * @brief Source of the high level state generator
 */

// Debug includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// Configuration includes
//
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/StateGenerator.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/SolveTiming/Before.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

   StateGenerator::StateGenerator()
      : SimulationBase()
   {
   }

   StateGenerator::~StateGenerator()
   {
   }

   void StateGenerator::preRun()
   {
   }

   void StateGenerator::mainRun()
   {
      // Solve trivial equations
      this->mPseudospectral.solveTrivialEquations(SolveTiming::Before::id());

      // Compute nonlinear terms
      this->mPseudospectral.computeNonlinear();

/// \mhdBug Problem with equations for generating exact states
      // Solve trivial equations
      this->mPseudospectral.solveTrivialEquations(SolveTiming::After::id());

      // Solve diagnostic equations
//      this->solveDiagnosticEquations(SolveTiming::After::id());

      // Synchronise computation nodes
      QuICCEnv().synchronize();
   }

   void StateGenerator::postRun()
   {
      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();

      // Write the output
      this->writeOutput();

      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();
   }

   void StateGenerator::tuneOutput()
   {
      // Get time and timestep from configuration
      Array tstep = this->mSimIoCtrl.config().timestepping();

      // Set to zero if value is negative
      if(tstep(0) < 0)
      {
         tstep(0) = 0;
      }
      if(tstep(1) < 0)
      {
         tstep(1) = 0;
      }

      // Loop over all files added to the simulation control
      SimulationIoControl::hdf5_iterator  fIt;
      for(fIt = this->mSimIoCtrl.beginHdf5(); fIt != this->mSimIoCtrl.endHdf5(); ++fIt)
      {
         (*fIt)->setSimTime(tstep(0), tstep(1));
      }
   }

   void StateGenerator::writeOutput()
   {
      // Write final state file (using stored time and timestep)
      this->mSimIoCtrl.writeHdf5(-1, -1);

      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();
   }

}
