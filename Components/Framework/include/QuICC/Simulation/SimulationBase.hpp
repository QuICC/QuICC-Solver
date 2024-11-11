/**
 * @file SimulationBase.hpp
 * @brief High level implementation of a base for the simulations
 */

#ifndef QUICC_SIMULATIONBASE_HPP
#define QUICC_SIMULATIONBASE_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Equations/EquationOptions.hpp"
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/Timers/ExecutionTimer.hpp"
#include "QuICC/SpatialScheme/IBuilder.hpp"
#include "QuICC/Simulation/SimulationRunControl.hpp"
#include "QuICC/Simulation/SimulationIoControl.hpp"
#include "QuICC/Simulation/SimulationBoundary.hpp"
#include "QuICC/Pseudospectral/Coordinator.hpp"
#include "QuICC/Io/Config/ConfigurationReader.hpp"
#include "QuICC/ScalarFields/ScalarField.hpp"
#include "QuICC/TypeSelectors/TransformCommSelector.hpp"
#include "QuICC/TypeSelectors/ParallelSelector.hpp"
#include "QuICC/TransformConfigurators/TransformTree.hpp"
#include "QuICC/LoadSplitter/LoadSplitter.hpp"
#include "QuICC/Io/Config/Simulation/Physical.hpp"
#include "QuICC/Io/Config/Simulation/Boundary.hpp"
#include "QuICC/Io/Variable/IVariableHdf5Reader.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"

namespace QuICC {

   /**
    * @brief High level implementation of a base for the simulations
    */
   class SimulationBase
   {
      public:
         /**
          * @brief Constructor
          */
         SimulationBase();

         /**
          * @brief Simple empty destructor
          */
         virtual ~SimulationBase() = default;

         /**
          * @brief Initialise the configuration read from file
          *
          * @param features   Features of model
          * @param modelHash  Version string of the model
          */
         void getConfig(std::map<std::string,MHDFloat>& cfg, std::set<SpatialScheme::Feature>& features, const std::string modelVersion);

         /**
          * @brief Initialise the configuration read from file
          */
         void updateConfig(const std::map<std::string,MHDFloat>& cfg);

         /**
          * @brief Initialise the base components of the simulation
          */
         void initBase();

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
          */
         template<typename TScheme> void initResolution(std::shared_ptr<TScheme> spScheme);

         /**
          * @brief Create the simulation wide boundary conditions
          */
         SharedSimulationBoundary createBoundary();

         /**
          * @brief Add equation to solver
          */
         template <typename TEquation> std::shared_ptr<TEquation> addEquation();

         /**
          * @brief Add equation to solver
          *
          * @param spOtions   Equation options
          */
         template <typename TEquation, typename TOptions> std::shared_ptr<TEquation> addEquation(std::shared_ptr<TOptions> spOptions);

         /**
          * @brief Add equation to solver
          *
          * @param spBackend Model backend
          */
         template <typename TEquation> std::shared_ptr<TEquation> addEquation(std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Add equation to solver
          *
          * @param spBackend Model backend
          * @param spBackend Equation options
          */
         template <typename TEquation, typename TOptions> std::shared_ptr<TEquation> addEquation(std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<TOptions> spOptions);

         /**
          * @brief Add graph description to solver
          * @param graphStr mlir module
          * @param physParams scaling parameters
          */
         void addGraph(const std::string& graphStr,
            const Graph::PhysicalParameters<MHDFloat>& physParams);

         /**
          * @brief Set the base simulation configuration file and parameters
          *
          * @param bcNames Vector of names for the boundary conditions
          */
         void setConfiguration(const int dimension, const std::string type, const std::vector<bool>& isPeriodicBox, const std::vector<std::string>& bcNames, const std::vector<std::string>& ndNames, const std::map<std::string,std::map<std::string,int> >& modelCfg);

         /**
          * @brief Set initial state through input file
          *
          * @param spInitFile Shared initial state file
          */
         void setInitialState(Io::Variable::SharedStateFileReader spInitFile);

         /**
          * @brief Add ASCII output file to solver
          *
          * @param spOutFile Shared ASCII output file
          */
         void addAsciiOutputFile(Io::Variable::SharedIVariableAsciiWriter spOutFile);

         /**
          * @brief Add HDF5 output file to solver
          *
          * @param spOutFile Shared HDF5 output file
          */
         void addHdf5OutputFile(Io::Variable::SharedIVariableHdf5NWriter spOutFile);

         /**
          * @brief Add Stats output file to solver
          *
          * @param spOutFile Shared stats output file
          */
         void addStatsOutputFile(Io::Stats::SharedIStatisticsAsciiWriter spOutFile);

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

         /**
          * @brief Get equation
          */
         const Equations::SharedEquationParameters eqParams() const;

      protected:
         /**
          * @brief Shared resolution
          */
         SharedResolution mspRes;

         /**
          * @brief Pseudospectral solver
          */
         Pseudospectral::Coordinator mPseudospectral;

         /**
          * @brief Shared resolution
          */
         Equations::SharedEquationParameters mspEqParams;

         /**
          * @brief Simulation run control
          */
         SimulationRunControl mSimRunCtrl;

         /**
          * @brief Simulation IO control
          */
         SimulationIoControl mSimIoCtrl;

      private:
         /**
          * @brief Add addition configuration file parts
          *
          * @param spCfgFile  Configuration file
          */
         virtual void addConfigurationNode(Io::Config::SharedConfigurationReader spCfgFile);

         /**
          * @brief Initialise the implementation specific base components
          */
         virtual void initAdditionalBase();

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

         /**
          * @brief Setup the output files added by the model
          */
         void setupOutput();

         /**
          * @brief Allow for implementation specific output tuning
          */
         virtual void tuneOutput();

         /**
          * @brief Allow for additional operators on the initial state input file
          *
          * @param spInitFile Shared initial state file
          */
         virtual void tuneInitialState(Io::Variable::SharedStateFileReader spInitFile);
   };

   template <typename TScheme> void SimulationBase::initResolution(std::shared_ptr<TScheme> spScheme)
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
      auto  spBuilder = spScheme->createBuilder(dim, true);
      splitter.init(spBuilder, {this->mSimIoCtrl.config().algorithm()}, this->mSimIoCtrl.config().grouper(), this->mSimIoCtrl.config().cpuFactors());

      stage.done();
      stage.start("extracting communication structure");

      // Get best splitting resolution object
      std::pair<SharedResolution, Parallel::SplittingDescription>  best = splitter.bestSplitting();

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

      // Initialize parallelisation setup for pseudospectral part
      this->mPseudospectral.initParallel(this->mspRes, best.second);
   }

   template <typename TEquation> std::shared_ptr<TEquation> SimulationBase::addEquation()
   {
      // Create a null backend
      auto spBackend = std::shared_ptr<Model::IModelBackend>();

      // Create default options
      auto spOptions = std::make_shared<Equations::EquationOptions>();

      return this->addEquation<TEquation>(spBackend);
   }

   template <typename TEquation, typename TOptions> std::shared_ptr<TEquation> SimulationBase::addEquation(std::shared_ptr<TOptions> spOptions)
   {
      // Create default backend
      auto spBackend = std::shared_ptr<Model::IModelBackend>();

      return this->addEquation<TEquation>(spBackend, spOptions);
   }

   template <typename TEquation> std::shared_ptr<TEquation> SimulationBase::addEquation(std::shared_ptr<Model::IModelBackend> spBackend)
   {
      // Create shared equation
      std::shared_ptr<TEquation> spEq = std::make_shared<TEquation>(this->mspEqParams, this->mspRes->sim().spSpatialScheme(), spBackend);

      // Create default options
      auto spOptions = std::make_shared<Equations::EquationOptions>();

      // Add shared equation
      this->mPseudospectral.addEquation(spEq, spOptions->it());

      return spEq;
   }

   template <typename TEquation, typename TOptions> std::shared_ptr<TEquation> SimulationBase::addEquation(std::shared_ptr<Model::IModelBackend> spBackend, std::shared_ptr<TOptions> spOptions)
   {
      // Create shared equation
      std::shared_ptr<TEquation> spEq = std::make_shared<TEquation>(this->mspEqParams, this->mspRes->sim().spSpatialScheme(), spBackend, spOptions);

      // Add shared equation
      this->mPseudospectral.addEquation(spEq, spOptions->it());

      return spEq;
   }

   /// Typedef for a shared pointer of a SimulationBase
   typedef std::shared_ptr<SimulationBase> SharedSimulationBase;

} // QuICC

#endif // QUICC_SIMULATIONBASE_HPP
