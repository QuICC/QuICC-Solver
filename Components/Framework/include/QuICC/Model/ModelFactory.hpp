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

      // Set spatial scheme
      auto spScheme = std::make_shared<typename TModel::SchemeType>(model.SchemeFormulation(), GridPurpose::SIMULATION);

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
      spSim->getConfig(cfg, features);
      spScheme->enable(features);
      model.configure(spScheme->features());
      spSim->updateConfig(model.backend().automaticParameters(cfg));

      spSim->initBase();


      // Initialise resolution
      spSim->initResolution(spScheme);

      // Add equations
      model.addEquations(spSim);

      // Add ASCII output files
      model.addAsciiOutputFiles(spSim);

      // Add HDF5 output files
      model.addHdf5OutputFiles(spSim);

      // Add statistics output files
      model.addStatsOutputFiles(spSim);

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
