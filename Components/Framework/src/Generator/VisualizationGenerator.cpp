/**
 * @file VisualizationGenerator.cpp
 * @brief Source of the high level state generator
 */

// First includes
//

// Debug includes
//

// Configuration includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/Profiler/ProfilerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Generator/VisualizationGenerator.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

   VisualizationGenerator::VisualizationGenerator()
      : SimulationBase()
   {
   }

   VisualizationGenerator::~VisualizationGenerator()
   {
   }

   void VisualizationGenerator::preRun()
   {
   }


   void VisualizationGenerator::mainRun()
   {
      // Solve the trivial equations
      this->mPseudospectral.explicitTrivialEquations(ModelOperator::ExplicitLinear::id());
      this->mPseudospectral.solveTrivialEquations(SolveTiming::After::id());

      // Compute nonlinear terms
      this->mPseudospectral.computeNonlinear();

      // Synchronise computation nodes
      QuICCEnv().synchronize();
   }

   void VisualizationGenerator::postRun()
   {
      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();

      // Write the output
      this->writeOutput();

      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();
   }

   void VisualizationGenerator::tuneOutput()
   {
   }

   void VisualizationGenerator::writeOutput()
   {
      // Write final state file (using store time and timestep)
      this->mSimIoCtrl.writeHdf5(this->mPseudospectral.startTime(), this->mPseudospectral.startTimestep());

      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();
   }

   void VisualizationGenerator::tuneInitialState(Io::Variable::SharedStateFileReader spInitFile)
   {
      // Get time from initial state
      MHDFloat time = spInitFile->time();

      // Get timestep from initial state
      MHDFloat timestep = spInitFile->timestep();

      // Loop over all files added to the simulation control
      SimulationIoControl::hdf5_iterator  fIt;
      for(fIt = this->mSimIoCtrl.beginHdf5(); fIt != this->mSimIoCtrl.endHdf5(); ++fIt)
      {
         (*fIt)->setSimTime(time, timestep);
      }
   }

}
