/**
 * @file ConfigFactory.hpp
 * @brief Implementation of the linear stability config factory
 */

#ifndef QUICC_CONFIGFACTORY_HPP
#define QUICC_CONFIGFACTORY_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/Io/Config/ConfigurationWriter.hpp"
#include "QuICC/Io/Config/Simulation/Boundary.hpp"
#include "QuICC/Io/Config/Simulation/Stability.hpp"
#include "QuICC/Io/Config/Simulation/Physical.hpp"
#include "QuICC/NonDimensional/GrowthRate.hpp"
#include "QuICC/NonDimensional/MaxIteration.hpp"
#include "QuICC/NonDimensional/Nev.hpp"
#include "QuICC/NonDimensional/Omega.hpp"
#include "QuICC/NonDimensional/Sort.hpp"
#include "QuICC/NonDimensional/StabilityMode.hpp"
#include "QuICC/NonDimensional/Tolerance.hpp"
#include "QuICC/NonDimensional/WriteMtx.hpp"
#include "QuICC/NonDimensional/registerStability.hpp"

namespace QuICC {

/**
 * @brief Implementation of the configuration for linear stability factory
 */
template <class TModel> class ConfigFactory
{
public:
   class ConfigApp
   {
   public:
      ConfigApp(std::shared_ptr<Io::Config::ConfigurationWriter> spWriter) :
          mspWriter(spWriter) {};

      /// Write configuration file
      void run()
      {
         this->mspWriter->init();
         this->mspWriter->write();
      };

      /// Finalize
      void finalize()
      {
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

   /**
    * @brief Constructor
    */
   ConfigFactory() = delete;

   /**
    * @brief Destructor
    */
   ~ConfigFactory() = delete;

protected:
   /**
    * @brief Add stability specific NonDimensional parameters
    */
   static void addParameters(std::vector<std::string>& ndNames);
};

template <class TModel>
void ConfigFactory<TModel>::addParameters(std::vector<std::string>& ndNames)
{
   // Add configuration parameters for Stability solver
   ndNames.push_back(NonDimensional::Tolerance().tag());
   ndNames.push_back(NonDimensional::MaxIteration().tag());
   ndNames.push_back(NonDimensional::Omega().tag());
   ndNames.push_back(NonDimensional::GrowthRate().tag());
   ndNames.push_back(NonDimensional::Nev().tag());
   ndNames.push_back(NonDimensional::Sort().tag());
   ndNames.push_back(NonDimensional::StabilityMode().tag());
   ndNames.push_back(NonDimensional::WriteMtx().tag());
}

template <class TModel>
typename ConfigFactory<TModel>::ReturnType ConfigFactory<TModel>::create()
{
   // Register stability specific nondimensional parameters
   NonDimensional::registerStability();

   // Create model
   TModel model;
   model.init();

   // Create and tune spatial scheme
   auto spScheme = std::make_shared<typename TModel::SchemeType>(
      model.SchemeFormulation(), QuICC::GridPurpose::SIMULATION);
   model.tuneScheme(spScheme);

   // Set dimension
   int dim = spScheme->dimension();

   // Set type string
   std::string type = spScheme->tag();

   // Box periodicity
   std::vector<bool> isPeriodicBox = model.backend().isPeriodicBox();

   // Create configuration writer
   auto writer = std::make_shared<QuICC::Io::Config::ConfigurationWriter>(dim,
      isPeriodicBox, type);

   // Create list of field ID strings for boundary conditions
   std::vector<std::string> bcNames = model.backend().fieldNames();

   // Create list of nondimensional ID strings for physical parameters
   std::vector<std::string> ndNames = model.backend().paramNames();
   ConfigFactory::addParameters(ndNames);

   // Add the physical part
   auto spPhys =
      std::make_shared<QuICC::Io::Config::Simulation::Physical>(ndNames);
   writer->rspSimulation()->addNode(QuICC::Io::Config::Simulation::PHYSICAL,
      spPhys);

   // Add the boundary part
   auto spBound =
      std::make_shared<QuICC::Io::Config::Simulation::Boundary>(bcNames);
   writer->rspSimulation()->addNode(QuICC::Io::Config::Simulation::BOUNDARY,
      spBound);

   // Add the stability field part
   std::vector<std::string> stabilityNames = {"stability_field", "stability_comp"};
   auto spStability =
      std::make_shared<QuICC::Io::Config::Simulation::Stability>(stabilityNames);
   writer->rspSimulation()->addNode(QuICC::Io::Config::Simulation::STABILITY,
      spStability);

   // Get model configuration tags
   auto modelCfg = model.configTags();

   // Add the model part
   writer->rspModel()->addNodes(modelCfg);

   auto spConfig = std::make_shared<ConfigFactory<TModel>::ConfigApp>(writer);

   return spConfig;
}

} // namespace QuICC

#endif // QUICC_CONFIGFACTORY_HPP
