/**
 * @file VisualizationGenerator.cpp
 * @brief Source of the high level state generator
 */

// Configuration includes
//
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// Project includes
//
#include "QuICC/Generator/VisualizationGenerator.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/SolveTiming/After.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Simulation/SimulationIoTools.hpp"

namespace QuICC {

   void VisualizationGenerator::preRun()
   {
      // Write ascii files
      SimulationIoTools::updateHeavyAscii(this->mSimIoCtrl.beginAscii(), this->mSimIoCtrl.endAscii(), this->mPseudospectral.transformCoordinator());
      this->mSimIoCtrl.writeAscii(this->mPseudospectral.startTime(), this->mPseudospectral.startTimestep());
   }

   void VisualizationGenerator::mainRun()
   {
      for(auto j: this->mPseudospectral.it())
      {
         // Solve the trivial equations
         this->mPseudospectral.explicitTrivialEquations(ModelOperator::ExplicitLinear::id(), j);
         this->mPseudospectral.solveTrivialEquations(SolveTiming::After::id(), j);

         // Compute nonlinear terms
         this->mPseudospectral.computeNonlinear(j);

         // Synchronise computation nodes
         QuICCEnv().synchronize();
      }
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

      // Write Hdf5 files
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
