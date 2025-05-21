/**
 * @file StabilityBase.cpp
 * @brief Source of the high level simulation base
 */

// System includes
//
#include <algorithm>
#include <stdexcept>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Bc/Scheme/Galerkin.hpp"
#include "QuICC/Debug/DebuggerMacro.h"
#include "QuICC/Debug/StorageProfiler/StorageProfilerMacro.h"
#include "QuICC/QuICCTimer.hpp"
#include "QuICC/Simulation/SimulationIoTools.hpp"
#include "QuICC/Tools/Formatter.hpp"
#include "Stability/StabilityBase.hpp"

namespace QuICC {

StabilityBase::StabilityBase()
{
   QuICCTimer().start();
}

SharedSimulationBoundary StabilityBase::createBoundary()
{
   // Create shared simulation boundary
   auto spBcs = std::make_shared<SimulationBoundary>(this->config().boundary());

   return spBcs;
}

void StabilityBase::setConfiguration(const int dimension,
   const std::string type, const std::vector<bool>& isPeriodicBox,
   const std::vector<std::string>& bcNames,
   const std::vector<std::string>& ndNames,
   const std::map<std::string, std::map<std::string, int>>& modelCfg)
{
   // Create shared configuration file
   auto spCfgFile = std::make_shared<Io::Config::ConfigurationReader>(dimension,
      isPeriodicBox, type);

   // Create the equation parameter shared pointer
   auto spEqParams = std::make_shared<Equations::EquationParameters>();
   this->mspEqParams = spEqParams;

   // Create the equation parameter dependent configuration part
   auto spPhys = std::make_shared<Io::Config::Simulation::Physical>(ndNames);

   // Add physical part to configuration file
   spCfgFile->rspSimulation()->addNode(Io::Config::Simulation::PHYSICAL,
      spPhys);

   // Create the boundary condition configuration part
   auto spBound = std::make_shared<Io::Config::Simulation::Boundary>(bcNames);

   // Add boundary part to configuration file
   spCfgFile->rspSimulation()->addNode(Io::Config::Simulation::BOUNDARY,
      spBound);

   // Add model part to configuration file
   spCfgFile->rspModel()->addNodes(modelCfg);

   // Add addition configuration parts
   this->addConfigurationNode(spCfgFile);

   // Set configuration file
   this->mSimIoCtrl.setConfigurationFile(spCfgFile);
}

void StabilityBase::getConfig(std::map<std::string, MHDFloat>& cfg,
   std::set<SpatialScheme::Feature>& features, const std::string modelVersion)
{
   // Initialise the IO system
   this->mSimIoCtrl.init(modelVersion);

   cfg = this->config().physical();

   if (this->config().bcScheme() == Bc::Scheme::Galerkin::id())
   {
      features.insert(SpatialScheme::Feature::GalerkinBasis);
   }
}

void StabilityBase::updateConfig(const std::map<std::string, MHDFloat>& cfg)
{
   this->mSimIoCtrl.rConfig().rPhysical().insert(cfg.begin(), cfg.end());
}

const Model::IModelBackend& StabilityBase::backend() const
{
   return *this->mspBackend;
}

void StabilityBase::initBase(std::shared_ptr<Model::IModelBackend> spBackend)
{
   this->mspBackend = spBackend;

   // Enable linearized equations
   this->mspBackend->enableLinearized(true);

   // Initialise the equation parameters
   this->mspEqParams->init(this->config().physical());

   // Get number CPU from configuration file
   int nCpu = this->config().nCpu();

   // Initialise the workflow
   QuICCEnv().setup(nCpu);

   // Make sure nodes are synchronised after initialisation
   QuICCEnv().synchronize();
}

void StabilityBase::init(const SharedSimulationBoundary spBcs)
{
   // Cleanup IO control
   this->mSimIoCtrl.cleanup();
}

void StabilityBase::run()
{
   // Stop timer and update initialisation time
   QuICCTimer().stop();
   QuICCTimer().update(ExecutionTimer::INIT);

   // Start timer
   QuICCTimer().start();

   // Execute pre-run steps
   this->preRun();

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

void StabilityBase::finalize()
{
   // Print execution timer infos
   QuICCTimer().printInfo(std::cout);

   // Print storage profiling infos (if required)
   StorageProfilerMacro_printInfo();
}

const SpatialScheme::ISpatialScheme& StabilityBase::ss() const
{
   return this->mspRes->sim().ss();
}

std::shared_ptr<const SpatialScheme::ISpatialScheme>
StabilityBase::spSpatialScheme() const
{
   return this->mspRes->sim().spSpatialScheme();
}

const SimulationConfig& StabilityBase::config() const
{
   return this->mSimIoCtrl.config();
}

void StabilityBase::addConfigurationNode(
   Io::Config::SharedConfigurationReader spCfgFile)
{
   // Empty to simplify implementations but can't be called from derived class
}
} // namespace QuICC
