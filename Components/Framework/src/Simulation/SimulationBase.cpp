/**
 * @file Simulation.cpp
 * @brief Source of the high level simulation base
 */

// Debug includes
//
#include "QuICC/Debug/DebuggerMacro.h"

// Configuration includes
//
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"

// System includes
//
#include <algorithm>
#include <stdexcept>

// Class include
//
#include "QuICC/Simulation/SimulationBase.hpp"

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/QuICCTimer.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "QuICC/Simulation/SimulationIoTools.hpp"
#include "QuICC/Bc/Scheme/Galerkin.hpp"

namespace QuICC {

   SimulationBase::SimulationBase()
      : mSimRunCtrl()
   {
      QuICCTimer().start();
   }

   SimulationBase::~SimulationBase()
   {
   }

   SharedSimulationBoundary SimulationBase::createBoundary()
   {
      // Create shared simulation boundary
      auto spBcs = std::make_shared<SimulationBoundary>(this->config().boundary());

      return spBcs;
   }

   void SimulationBase::setConfiguration(const int dimension, const std::string type, const std::vector<bool>& isPeriodicBox, const std::vector<std::string>& bcNames, const std::vector<std::string>& ndNames, const std::map<std::string,std::map<std::string,int> >& modelCfg)
   {
      // Create shared configuration file
      auto spCfgFile = std::make_shared<Io::Config::ConfigurationReader>(dimension, isPeriodicBox, type);

      // Create the equation parameter shared pointer
      auto spEqParams = std::make_shared<Equations::EquationParameters>();
      this->mspEqParams = spEqParams;

      // Create the equation parameter dependent configuration part
      auto spPhys = std::make_shared<Io::Config::Simulation::Physical>(ndNames);

      // Add physical part to configuration file
      spCfgFile->rspSimulation()->addNode(Io::Config::Simulation::PHYSICAL, spPhys);

      // Create the boundary condition configuration part
      auto spBound = std::make_shared<Io::Config::Simulation::Boundary>(bcNames);

      // Add boundary part to configuration file
      spCfgFile->rspSimulation()->addNode(Io::Config::Simulation::BOUNDARY, spBound);

      // Add model part to configuration file
      spCfgFile->rspModel()->addNodes(modelCfg);

      // Add addition configuration parts
      this->addConfigurationNode(spCfgFile);

      // Set configuration file
      this->mSimIoCtrl.setConfigurationFile(spCfgFile);
   }

   void SimulationBase::getConfig(std::map<std::string,MHDFloat>& cfg, std::set<SpatialScheme::Feature>& features, const std::string modelVersion)
   {
      // Initialise the IO system
      this->mSimIoCtrl.init(modelVersion);

      cfg = this->config().physical();

      if(this->config().bcScheme() == Bc::Scheme::Galerkin::id())
      {
         features.insert(SpatialScheme::Feature::GalerkinBasis);
      }
   }

   void SimulationBase::updateConfig(const std::map<std::string,MHDFloat>& cfg)
   {
      this->mSimIoCtrl.rConfig().rPhysical().insert(cfg.begin(), cfg.end());
   }

   void SimulationBase::initBase()
   {
      // Initialise the equation parameters
      this->mspEqParams->init(this->config().physical());

      // Get number CPU from configuration file
      int nCpu = this->config().nCpu();

      // Initialise the workflow
      QuICCEnv().setup(nCpu);

      // Initialise additional things depending on implementation
      this->initAdditionalBase();

      // Make sure nodes are synchronised after initialisation
      QuICCEnv().synchronize();
   }

   void SimulationBase::init(const SharedSimulationBoundary spBcs)
   {
      // Get timestep information from configuration file
      Array tstep = this->config().timestepping();

      this->mPseudospectral.init(tstep, spBcs);

      StageTimer stage;
      stage.start("setup output files");

      // Setup output files (ASCII diagnostics, state files, etc)
      this->setupOutput();

      stage.done();

      // Cleanup IO control
      this->mSimIoCtrl.cleanup();
   }

   void SimulationBase::run()
   {
      // Final initialisation of the solvers
      this->mPseudospectral.initSolvers();

      // Stop timer and update initialisation time
      QuICCTimer().stop();
      QuICCTimer().update(ExecutionTimer::INIT);

      // Start timer
      QuICCTimer().start();

      // Execute pre-run steps
      this->preRun();

      // Final initialisation of the solvers
      this->mPseudospectral.cleanupForRun();

      // Stop pre-run timing
      QuICCTimer().stop();
      QuICCTimer().update(ExecutionTimer::PRERUN);

      // Synchronize computation nodes
      QuICCEnv().synchronize();
      StageTimer::completed("Simulation initialisation successfull");

      // Start timer
      QuICCTimer().start();

      // Do main loop
      this->mainRun();

      // Stop main loop timing
      QuICCTimer().stop();
      QuICCTimer().update(ExecutionTimer::RUN);

      // Start timer
      QuICCTimer().start();

      // Execute post-run operations
      this->postRun();

      // Stop post-run timing
      QuICCTimer().stop();
      QuICCTimer().update(ExecutionTimer::POSTRUN);

      // Synchronise computation nodes
      QuICCEnv().synchronize();
   }

   void SimulationBase::finalize()
   {
      // Print simulation run infos
      this->mSimRunCtrl.printInfo(std::cout);

      // Print execution timer infos
      QuICCTimer().printInfo(std::cout);

      // Print storage profiling infos (if required)
      StorageProfilerMacro_printInfo();
   }

   void SimulationBase::setInitialState(Io::Variable::SharedStateFileReader spInitFile)
   {
      // Loop over all scalars
      for(auto scalIt = this->mPseudospectral.scalarVariables().begin(); scalIt != this->mPseudospectral.scalarVariables().end(); scalIt++)
      {
         spInitFile->addScalar((*scalIt));
      }

      // Loop over all vector variables
      for(auto vectIt = this->mPseudospectral.vectorVariables().begin(); vectIt != this->mPseudospectral.vectorVariables().end(); vectIt++)
      {
         spInitFile->addVector((*vectIt));
      }

      // Loop over all imposed scalars
      for(auto scalIt = this->mPseudospectral.imposedScalarVariables().begin(); scalIt != this->mPseudospectral.imposedScalarVariables().end(); scalIt++)
      {
         spInitFile->addScalar((*scalIt));
      }

      // Loop over all imposed vector variables
      for(auto vectIt = this->mPseudospectral.imposedVectorVariables().begin(); vectIt != this->mPseudospectral.imposedVectorVariables().end(); vectIt++)
      {
         spInitFile->addVector((*vectIt));
      }

      // Initialise file
      spInitFile->init();

      // Read in data
      spInitFile->read();

      // Addition operations on initial state file
      this->tuneInitialState(spInitFile);

      // Forward state file time and timestep to diagnostic coordinator
      this->mPseudospectral.useStateTime(spInitFile->time(), spInitFile->timestep());

      // Finalise file
      spInitFile->finalize();
   }

   void SimulationBase::tuneInitialState(Io::Variable::SharedStateFileReader spInitFile)
   {
   }

   void SimulationBase::addAsciiOutputFile(Io::Variable::SharedIVariableAsciiWriter spOutFile)
   {
      this->mSimIoCtrl.addAsciiOutputFile(spOutFile);
   }

   void SimulationBase::addHdf5OutputFile(Io::Variable::SharedIVariableHdf5NWriter spOutFile)
   {
      this->mSimIoCtrl.addHdf5OutputFile(spOutFile);
   }

   void SimulationBase::addStatsOutputFile(Io::Stats::SharedIStatisticsAsciiWriter spOutFile)
   {
      this->mSimIoCtrl.addStatsOutputFile(spOutFile);
   }

   void SimulationBase::setupOutput()
   {
      // Loop over all ASCII files added to the simulation control
      SimulationIoControl::ascii_iterator  asciiIt;
      for(asciiIt = this->mSimIoCtrl.beginAscii(); asciiIt != this->mSimIoCtrl.endAscii(); ++asciiIt)
      {
         // Loop over all scalars
         for(auto scalIt = this->mPseudospectral.scalarVariables().begin(); scalIt != this->mPseudospectral.scalarVariables().end(); scalIt++)
         {
            (*asciiIt)->addScalar((*scalIt));
         }

         // Loop over all vector variables
         for(auto vectIt = this->mPseudospectral.vectorVariables().begin(); vectIt != this->mPseudospectral.vectorVariables().end(); vectIt++)
         {
            (*asciiIt)->addVector((*vectIt));
         }

         // Set mesh
         if((*asciiIt)->space() == Dimensions::Space::PHYSICAL)
         {
            (*asciiIt)->setMesh(this->mPseudospectral.transformCoordinator().mesh());
         }
      }

      // Loop over all HDF5 files added to the simulation control
      SimulationIoControl::hdf5_iterator  hdf5It;
      for(hdf5It = this->mSimIoCtrl.beginHdf5(); hdf5It != this->mSimIoCtrl.endHdf5(); ++hdf5It)
      {
         // Loop over all scalars
         for(auto scalIt = this->mPseudospectral.scalarVariables().begin(); scalIt != this->mPseudospectral.scalarVariables().end(); scalIt++)
         {
            (*hdf5It)->addScalar((*scalIt));
         }

         // Loop over all vector variables
         for(auto vectIt = this->mPseudospectral.vectorVariables().begin(); vectIt != this->mPseudospectral.vectorVariables().end(); vectIt++)
         {
            (*hdf5It)->addVector((*vectIt));
         }

         // Loop over all imposed scalars
         for(auto scalIt = this->mPseudospectral.imposedScalarVariables().begin(); scalIt != this->mPseudospectral.imposedScalarVariables().end(); scalIt++)
         {
            (*hdf5It)->addScalar((*scalIt));
         }

         // Loop over all imposed vector variables
         for(auto vectIt = this->mPseudospectral.imposedVectorVariables().begin(); vectIt != this->mPseudospectral.imposedVectorVariables().end(); vectIt++)
         {
            (*hdf5It)->addVector((*vectIt));
         }

         // Set physical mesh
         if((*hdf5It)->space() == Dimensions::Space::PHYSICAL)
         {
            (*hdf5It)->setMesh(this->mPseudospectral.transformCoordinator().mesh());
         }
      }

      // Loop over all satistics files added to the simulation control
      SimulationIoControl::stats_iterator  statsIt;
      for(statsIt = this->mSimIoCtrl.beginStats(); statsIt != this->mSimIoCtrl.endStats(); ++statsIt)
      {
         // Loop over all scalars
         for(auto scalIt = this->mPseudospectral.scalarVariables().begin(); scalIt != this->mPseudospectral.scalarVariables().end(); scalIt++)
         {
            (*statsIt)->addScalar((*scalIt));
         }

         // Loop over all vector variables
         for(auto vectIt = this->mPseudospectral.vectorVariables().begin(); vectIt != this->mPseudospectral.vectorVariables().end(); vectIt++)
         {
            (*statsIt)->addVector((*vectIt));
         }

         (*statsIt)->setMesh(this->mPseudospectral.transformCoordinator().mesh());
      }

      // Allow for implementation specific tuning
      this->tuneOutput();

      // init the output writers
      this->mSimIoCtrl.initWriters();
   }

   const SpatialScheme::ISpatialScheme& SimulationBase::ss() const
   {
      return this->mspRes->sim().ss();
   }

   std::shared_ptr<const SpatialScheme::ISpatialScheme> SimulationBase::spSpatialScheme() const
   {
      return this->mspRes->sim().spSpatialScheme();
   }

   const SimulationConfig& SimulationBase::config() const
   {
      return this->mSimIoCtrl.config();
   }

   void SimulationBase::tuneOutput()
   {
   }

   void SimulationBase::addConfigurationNode(Io::Config::SharedConfigurationReader spCfgFile)
   {
      // Empty to simplify implementations but can't be called from derived class
   }

   void SimulationBase::initAdditionalBase()
   {
      // Empty to simplify implementations but can't be called from derived class
   }
}
