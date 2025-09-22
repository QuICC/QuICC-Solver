/** 
 * @file StateGeneratorFactory.hpp
 * @brief Implementation of the state generator model factory
 */

#ifndef QUICC_STATEGENERATORFACTORY_HPP
#define QUICC_STATEGENERATORFACTORY_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/Model/IPhysicalModel.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the state generator model factory 
    */
   template <class TModel> class StateGeneratorFactory
   {
      public:
         /// Typdef of type of oject created
         typedef std::shared_ptr<StateGenerator> ReturnType;

         /**
          * @brief Create a shared state generator for the model
          */
         static ReturnType create();

         /**
          * @brief Constructor
          */
         StateGeneratorFactory() = delete;

         /**
          * @brief Destructor
          */
         ~StateGeneratorFactory() = delete;
      protected:

      private:
   };

   template <class TModel> typename StateGeneratorFactory<TModel>::ReturnType StateGeneratorFactory<TModel>::create()
   {
      // Create model
      TModel model;
      model.init();

      // Create simulation
      auto spGen = std::make_shared<StateGenerator>();

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
      spGen->setConfiguration(spScheme->dimension(), spScheme->tag(), model.backend().isPeriodicBox(), bcNames, ndNames, modelCfg);

      // Initialise simulation
      std::map<std::string,MHDFloat> cfg;
      std::set<SpatialScheme::Feature> features;
      spGen->getConfig(cfg, features, model.version());
      spScheme->enable(features);
      model.configure(spScheme->features());
      spGen->updateConfig(model.backend().automaticParameters(cfg));

      spGen->initBase();

      // Initialise resolution
      spGen->initResolution(spScheme);

      // Add initial state generation equations
      model.addStates(spGen);

      // Set the boundary conditions
      SharedSimulationBoundary spBcs = spGen->createBoundary();

      // Initialise the simulation
      spGen->init(spBcs);

      // Set initial state
      model.setGeneratorState(spGen);

      return spGen;
   }

}

#endif // QUICC_STATEGENERATORFACTORY_HPP
