/** 
 * @file ConfigFactory.hpp
 * @brief Implementation of the physical model configuration factory
 */

#ifndef QUICC_CONFIGFACTORY_HPP
#define QUICC_CONFIGFACTORY_HPP

// System includes
//
#include <memory>
#include <vector>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Config/ConfigurationWriter.hpp"
#include "QuICC/Io/Config/Simulation/Physical.hpp"
#include "QuICC/Io/Config/Simulation/Boundary.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the physical model factory 
    */
   template <class TModel> class ConfigFactory
   {
      public:
         class ConfigApp
         {
            public:
               ConfigApp(std::shared_ptr<Io::Config::ConfigurationWriter> spWriter) 
                  : mspWriter(spWriter)
               {};

               /// Write configuration file
               void run() {
                  this->mspWriter->init();
                  this->mspWriter->write();
               };

               /// Finalize
               void finalize() {
                  this->mspWriter->finalize();
               };

            private:
               std::shared_ptr<Io::Config::ConfigurationWriter> mspWriter;
         };

         /// Typedef of type of object created
         typedef std::shared_ptr<ConfigApp> ReturnType;

         /**
          * @brief Create a shared simulation for the model
          */
         static ReturnType create();

      protected:

      private:
         /**
          * @brief Constructor
          */
         ConfigFactory();

         /**
          * @brief Destructor
          */
         ~ConfigFactory();
   };

   template <class TModel> typename ConfigFactory<TModel>::ReturnType ConfigFactory<TModel>::create()
   {
      // create model
      TModel model;
      model.init();

      // Create and tune spatial scheme
      auto spScheme = std::make_shared<typename TModel::SchemeType>(model.SchemeFormulation(), GridPurpose::SIMULATION);
      model.tuneScheme(spScheme);

      // Set dimension
      int dim = spScheme->dimension();

      // Set type string
      std::string type = spScheme->tag();

      // Box periodicity
      std::vector<bool> isPeriodicBox = model.backend().isPeriodicBox();

      // Create configuration writer
      auto writer = std::make_shared<Io::Config::ConfigurationWriter>(dim, isPeriodicBox, type);

      // Create list of field ID strings for boundary conditions
      std::vector<std::string> bcNames = model.backend().fieldNames();

      // Create list of nondimensional ID strings for physical parameters
      std::vector<std::string> ndNames = model.backend().paramNames();

      // Add the physical part
      auto spPhys = std::make_shared<Io::Config::Simulation::Physical>(ndNames);
      writer->rspSimulation()->addNode(Io::Config::Simulation::PHYSICAL, spPhys);

      // Add the boundary part
      auto spBound = std::make_shared<Io::Config::Simulation::Boundary>(bcNames);
      writer->rspSimulation()->addNode(Io::Config::Simulation::BOUNDARY, spBound);

      // Get model configuration tags
      auto modelCfg = model.configTags();

      // Add the model part
      writer->rspModel()->addNodes(modelCfg);

      auto spConfig = std::make_shared<ConfigFactory<TModel>::ConfigApp>(writer);

      return spConfig;
   }

}

#endif // QUICC_CONFIGFACTORY_HPP
