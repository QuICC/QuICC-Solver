/**
 * @file Simulation.cpp
 * @brief Source of the high level simulation
 */

// First includes
//
#include "QuICC/Equations/Tools/EquationTools.hpp"

// System includes
//
#include <algorithm>

// Project includes
//
#include "QuICC/Simulation/Simulation.hpp"
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/QuICCTimer.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/RuntimeStatus/GoOn.hpp"
#include "QuICC/Simulation/SimulationIoTools.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "Profiler/Interface.hpp"

namespace QuICC {

   Simulation::Simulation()
      : SimulationBase()
   {
   }

   Simulation::~Simulation()
   {
   }

   void Simulation::initAdditionalBase()
   {
      // Get the run configuration
      Array cfgRun = this->mSimIoCtrl.config().run();

      // Set the maximum simulation time
      this->mSimRunCtrl.setMaxSimTime(cfgRun(0));

      // Set the maximum wall time
      this->mSimRunCtrl.setMaxWallTime(cfgRun(1));
   }

   void Simulation::mainRun()
   {
      Profiler::RegionFixture<1> simFix("Simulation::mainRun");

      StageTimer::stage("Starting simulation");

      // Start main loop of simulation
      while(this->mSimRunCtrl.status() == RuntimeStatus::GoOn::id())
      {
         // Evolve pseudospectral equations
         this->mPseudospectral.evolve();

         // Update simulation run control
         this->mSimRunCtrl.updateSimulation(this->mPseudospectral.time(), this->mPseudospectral.timestep());

         // Update simulation IO control
         this->mSimIoCtrl.update();

         // Write the output
         this->writeOutput();

         // Synchronise computation nodes
         QuICCEnv().synchronize();

         // Update simulation run control
         this->mSimRunCtrl.updateCluster(QuICCTimer().queryTime(ExecutionTimer::TOTAL));

         // Increase iteration counter for timer
         QuICCTimer().iteration();
      }

      // Profile storage
      this->mPseudospectral.profileStorage();
   }

   void Simulation::preRun()
   {
      Profiler::RegionFixture<1> simFix("Simulation::preRun");

      StageTimer stage;
      stage.start("Prepare time evolution");

      this->mPseudospectral.prepareEvolution(this->mSimIoCtrl.config().timestepper());

      // Print timestepper information
      this->mPseudospectral.printInfo(std::cout);

      stage.done();
      stage.start("write initial ASCII files");

      // Update heavy calculation required for ASCII output
      SimulationIoTools::updateHeavyAscii(this->mSimIoCtrl.beginAscii(), this->mSimIoCtrl.endAscii(), this->mPseudospectral.transformCoordinator());

      // Write initial ASCII output
      this->mSimIoCtrl.writeAscii(this->mPseudospectral.time(), this->mPseudospectral.timestep());

      stage.done();
      stage.start("write initial HDF5 files");

      // Write initial state file
      this->mSimIoCtrl.writeHdf5(this->mPseudospectral.time(), this->mPseudospectral.timestep());

      stage.done();
      stage.start("write initial statistics files");

      // Update calculation required for statistics output
      this->mSimIoCtrl.prepareStats(this->mPseudospectral.time(), this->mPseudospectral.timestep());
      SimulationIoTools::updateStatsPre(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mPseudospectral.transformCoordinator());
      SimulationIoTools::updateStats(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mPseudospectral.transformCoordinator());
      SimulationIoTools::updateStatsPost(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mPseudospectral.transformCoordinator());

      // Write initial statistics output
      this->mSimIoCtrl.writeStats();

      stage.done();
   }

   void Simulation::writeOutput()
   {
      Profiler::RegionFixture<1> fix("Simulation::writeOutput");

      this->mSimIoCtrl.writeStats();

      if(this->mSimIoCtrl.isAsciiTime())
      {
         // Update heavy calculation required for ASCII output
         SimulationIoTools::updateHeavyAscii(this->mSimIoCtrl.beginAscii(), this->mSimIoCtrl.endAscii(), this->mPseudospectral.transformCoordinator());
      }

      // Write initial ASCII and HDF5 output files if applicable
      this->mSimIoCtrl.writeFiles(this->mPseudospectral.time(), this->mPseudospectral.timestep());

      // Write equation diagnotics
      this->mPseudospectral.writeDiagnostics(this->mSimIoCtrl.isAsciiTime(), this->mSimIoCtrl.isHdf5Time());
   }

   void Simulation::postRun()
   {
      Profiler::RegionFixture<1> simFix("Simulation::postRun");

      this->mSimIoCtrl.writeStats();

      StageTimer::stage("Post simulation");
      StageTimer  stage;

      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();

      stage.start("write final ASCII files");

      // Update heavy calculation required for ASCII output
      SimulationIoTools::updateHeavyAscii(this->mSimIoCtrl.beginAscii(), this->mSimIoCtrl.endAscii(), this->mPseudospectral.transformCoordinator());

      // Write final ASCII output
      this->mSimIoCtrl.writeAscii(this->mPseudospectral.time(), this->mPseudospectral.timestep());

      stage.done();
      stage.start("write final HDF5 files");

      // Write final state file
      this->mSimIoCtrl.writeHdf5(this->mPseudospectral.time(), this->mPseudospectral.timestep());

      stage.done();
      stage.start("write final statistics files");

//      // Update calculation required for statistics output
//      this->mSimIoCtrl.prepareStats(this->mPseudospectral.time(), this->mPseudospectral.timestep());
//      SimulationIoTools::updateStatsPre(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mPseudospectral.transformCoordinator());
//      SimulationIoTools::updateStats(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mPseudospectral.transformCoordinator());
//      SimulationIoTools::updateStatsPost(this->mSimIoCtrl.beginStats(), this->mSimIoCtrl.endStats(), this->mPseudospectral.transformCoordinator());
//
//      // Write final statistics output
//      this->mSimIoCtrl.writeStats();

      stage.done();

      // Synchronise all nodes of simulation
      QuICCEnv().synchronize();
   }

}
