/**
 * @file StabilityBase.hpp
 * @brief High level implementation of a stability solver base
 */

#ifndef QUICC_STABILITYBASE_HPP
#define QUICC_STABILITYBASE_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Equations/EquationParameters.hpp"
#include "QuICC/Io/Config/ConfigurationReader.hpp"
#include "QuICC/Io/Config/Simulation/Boundary.hpp"
#include "QuICC/Io/Config/Simulation/Physical.hpp"
#include "QuICC/LoadSplitter/LoadSplitter.hpp"
#include "QuICC/Model/IModelBackend.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/Simulation/SimulationBoundary.hpp"
#include "QuICC/Simulation/SimulationIoControl.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"
#include "QuICC/Timers/ExecutionTimer.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/TypeSelectors/ParallelSelector.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"

namespace QuICC {

/**
 * @brief High level implementation of a stability solver base
 */
class StabilityBase
{
public:
   /**
    * @brief Constructor
    */
   StabilityBase();

   /**
    * @brief Simple empty destructor
    */
   virtual ~StabilityBase() = default;

   /**
    * @brief Initialise the configuration read from file
    *
    * @param cfg           Configuration map
    * @param features      Features of model
    * @param modelVersion  Version string of the model
    */
   void getConfig(std::map<std::string, MHDFloat>& cfg,
      std::set<SpatialScheme::Feature>& features,
      const std::string modelVersion);

   /**
    * @brief Initialise the configuration read from file
    *
    * @param cfg           Configuration map
    */
   void updateConfig(const std::map<std::string, MHDFloat>& cfg);

   /**
    * @brief Initialise the base components of the simulation
    *
    * @param spBackend   Model backend
    */
   void initBase(std::shared_ptr<Model::IModelBackend> spBackend);

   /**
    * @brief Initialise the different components of the simulation
    *
    * @param spBcs Boundary condition information
    */
   void init(const SharedSimulationBoundary spBcs);

   /**
    * @brief Run the simulation
    */
   void run();

   /**
    * @brief Finalise simulation run
    */
   void finalize();

   /**
    * @brief Initialise the resolution
    *
    * @param spScheme   Spatial scheme
    */
   template <typename TScheme>
   void initResolution(std::shared_ptr<TScheme> spScheme);

   /**
    * @brief Create the simulation wide boundary conditions
    */
   SharedSimulationBoundary createBoundary();

   /**
    * @brief Set the base simulation configuration file and parameters
    *
    * @param dimension     Number of dimensions
    * @param type          Spatial scheme type
    * @param isPeriodicBox    Geometry is periodic?
    * @param bcNames Vector of names for the boundary conditions
    * @param ndNames Vector of names of nondimensional parameters
    * @param modelCfg   Model configuration
    */
   void setConfiguration(const int dimension, const std::string type,
      const std::vector<bool>& isPeriodicBox,
      const std::vector<std::string>& bcNames,
      const std::vector<std::string>& ndNames,
      const std::map<std::string, std::map<std::string, int>>& modelCfg);

   /**
    * @brief Get SpatialScheme
    */
   const SpatialScheme::ISpatialScheme& ss() const;

   /**
    * @brief Get SpatialScheme
    */
   std::shared_ptr<const SpatialScheme::ISpatialScheme> spSpatialScheme() const;

   /**
    * @brief Get configuration
    */
   const SimulationConfig& config() const;

protected:
   /**
    * @brief Get model backend
    */
   const Model::IModelBackend& backend() const;

   /**
    * @brief Shared resolution
    */
   SharedResolution mspRes;

   /**
    * @brief Shared resolution
    */
   Equations::SharedEquationParameters mspEqParams;

   /**
    * @brief Simulation IO control
    */
   SimulationIoControl mSimIoCtrl;

   /**
    * @brief model backend
    */
   std::shared_ptr<Model::IModelBackend> mspBackend;

private:
   /**
    * @brief Add addition configuration file parts
    *
    * @param spCfgFile  Configuration file
    */
   virtual void addConfigurationNode(
      Io::Config::SharedConfigurationReader spCfgFile);

   /**
    * @brief Do operations required just before starting the main loop
    */
   virtual void preRun() = 0;

   /**
    * @brief Do operations required during the main loop
    */
   virtual void mainRun() = 0;

   /**
    * @brief Do operations required just after finishing the main loop
    */
   virtual void postRun() = 0;
};

template <typename TScheme>
void StabilityBase::initResolution(std::shared_ptr<TScheme> spScheme)
{
   StageTimer stage;
   stage.start("optimizing load distribution");

   // Create the load splitter
   Parallel::LoadSplitter splitter(QuICCEnv().id(), QuICCEnv().size());

   // Extract dimensions from configuration file
   ArrayI dim = this->mSimIoCtrl.config().dimension();

   // Select transform implementations
   spScheme->setImplementation(this->mSimIoCtrl.config().transformSetup());

   // Initialise the load splitter
   auto spBuilder = spScheme->createBuilder(dim, true);
   splitter.init(spBuilder, {this->mSimIoCtrl.config().algorithm()},
      this->mSimIoCtrl.config().grouper());

   stage.done();
   stage.start("extracting communication structure");

   // Get best splitting resolution object
   std::pair<SharedResolution, Parallel::SplittingDescription> best =
      splitter.bestSplitting();

   // Store the shared resolution object
   this->mspRes = best.first;

   // Set additional options on final resolution object
   spBuilder->tuneResolution(this->mspRes, best.second);

   stage.done();

   // Extract box scale from configuration file
   Array box = this->mSimIoCtrl.config().boxScale();

   // Set the box scale
   this->mspRes->setBoxScale(box);
   // Set spatial scheme
   this->mspRes->setSpatialScheme(spScheme);
}

/// Typedef for a shared pointer of a StabilityBase
typedef std::shared_ptr<StabilityBase> SharedStabilityBase;

} // namespace QuICC

#endif // QUICC_STABILITYBASE_HPP
