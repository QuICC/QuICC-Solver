/** 
 * @file VisualizationGeneratorFactory.hpp
 * @brief Implementation of the visualization generator model factory
 */

#ifndef QUICC_VISUALIZATIONGENERATORFACTORY_HPP
#define QUICC_VISUALIZATIONGENERATORFACTORY_HPP

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
#include "QuICC/Generator/VisualizationGenerator.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the visualization generator model factory 
    */
   template <class TModel> class VisualizationGeneratorFactory
   {
      public:
         /**
          * @brief Create a shared state generator for the model
          */
         static SharedVisualizationGenerator createVisualization();

      protected:

      private:
         /**
          * @brief Constructor
          */
         VisualizationGeneratorFactory();

         /**
          * @brief Destructor
          */
         ~VisualizationGeneratorFactory();
   };

   template <class TModel> SharedVisualizationGenerator VisualizationGeneratorFactory<TModel>::createVisualization()
   {
      // Create model
      TModel model;
      model.init();

      // Create simulation
      auto spVis = std::make_shared<VisualizationGenerator>();

      // Set spatial scheme
      auto spScheme = std::make_shared<typename TModel::SchemeType>(model.SchemeFormulation(), GridPurpose::VISUALIZATION);

      // Create list of field ID strings for boundary conditions
      std::vector<std::string> bcNames = model.backend().fieldNames();

      // Create list of nondimensional ID strings for physical parameters
      std::vector<std::string> ndNames = model.backend().paramNames();

      // Get model configuration tags
      auto modelCfg = model.configTags();

      // Add configuration file and parameters
      spVis->setConfiguration(spScheme->dimension(), spScheme->tag(), model.backend().isPeriodicBox(), bcNames, ndNames, modelCfg);

      // Initialise simulation
      std::map<std::string,MHDFloat> cfg;
      std::set<SpatialScheme::Feature> features;
      spVis->getConfig(cfg, features);
      spScheme->enable(features);
      model.configure(spScheme->features());
      spVis->updateConfig(model.backend().automaticParameters(cfg));

      spVis->initBase();

      // Initialise resolution
      spVis->initResolution(spScheme);

      // Add initial state generation equations
      model.addVisualizers(spVis);

      // Set the boundary conditions
      SharedSimulationBoundary spBcs = spVis->createBoundary();

      // Initialise the simulation
      spVis->init(spBcs);

      // Set initial state
      model.setVisualizationState(spVis);

      return spVis;
   }

}

#endif // QUICC_VISUALIZATIONGENERATORFACTORY_HPP
