/** 
 * @file ModelFactory.hpp
 * @brief Implementation of the physical model factory
 */

#ifndef QUICC_MODELFACTORY_HPP
#define QUICC_MODELFACTORY_HPP

// First include
//

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/Simulation/Simulation.hpp"
#include "QuICC/Model/IPhysicalModel.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the physical model factory 
    */
   template <class TModel> class ModelFactory
   {
      public:
         /**
          * @brief Create a shared simulation for the model
          */
         static SharedSimulation createSimulation();

      protected:

      private:
         /**
          * @brief Constructor
          */
         ModelFactory();

         /**
          * @brief Destructor
          */
         ~ModelFactory();
   };

   template <class TModel> SharedSimulation ModelFactory<TModel>::createSimulation()
   {
      // Create model
      TModel model;
      model.init();

      // Create simulation
      auto spSim = std::make_shared<Simulation>();

      // Set and tune spatial scheme
      auto spScheme = std::make_shared<typename TModel::SchemeType>(model.SchemeFormulation(), GridPurpose::SIMULATION);
      model.tuneScheme(spScheme);

      // Create list of field ID strings for boundary conditions
      std::vector<std::string> bcNames = model.backend().fieldNames();

      // Create list of nondimensional ID strings for physical parameters
      std::vector<std::string> ndNames = model.backend().paramNames();

      // Get model configuration tags
      auto modelCfg = model.configTags();

      // Add configuration file and parameters
      spSim->setConfiguration(spScheme->dimension(), spScheme->tag(), model.backend().isPeriodicBox(), bcNames, ndNames, modelCfg);

      // Initialise simulation
      std::map<std::string,MHDFloat> cfg;
      std::set<SpatialScheme::Feature> features;
      spSim->getConfig(cfg, features, model.version());
      spScheme->enable(features);
      model.configure(spScheme->features());
      spSim->updateConfig(model.backend().automaticParameters(cfg));

      spSim->initBase();


      // Initialise resolution
      spSim->initResolution(spScheme);

      StageTimer stage;

      // Add equations
      stage.start("adding model equations");
      model.addEquations(spSim);
      stage.done();

      // Add ASCII output files
      stage.start("adding model ASCII output");
      model.addAsciiOutputFiles(spSim);
      stage.done();

      // Add HDF5 output files
      stage.start("adding model HDF5 output");
      model.addHdf5OutputFiles(spSim);
      stage.done();

      // Add statistics output files
      stage.start("adding model Stats output");
      model.addStatsOutputFiles(spSim);
      stage.done();

      // Set the boundary conditions
      SharedSimulationBoundary spBcs = spSim->createBoundary();

      // Initialise the simulation
      spSim->init(spBcs);

      // Set initial state
      model.setInitialState(spSim);

      return spSim;
   }

}

#endif // QUICC_MODELFACTORY_HPP
