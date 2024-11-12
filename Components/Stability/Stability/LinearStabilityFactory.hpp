/**
 * @file LinearStabilityFactory.hpp
 * @brief Implementation of the linear stability factory
 */

#ifndef QUICC_LINEARSTABILITYFACTORY_HPP
#define QUICC_LINEARSTABILITYFACTORY_HPP

// System includes
//

// Project includes
//
#include "QuICC/Enums/GridPurpose.hpp"
#include "QuICC/Model/IPhysicalModel.hpp"
#include "Stability/MarginalCurve.hpp"

namespace QuICC {

/**
 * @brief Implementation of the linear stability factory
 */
template <class TModel> class LinearStabilityFactory
{
public:
   /**
    * @brief Create a shared simulation for the model
    */
   static SharedMarginalCurve createSolver();

protected:
private:
   /**
    * @brief Constructor
    */
   LinearStabilityFactory() = default;

   /**
    * @brief Destructor
    */
   ~LinearStabilityFactory() = default;
};

template <class TModel>
SharedMarginalCurve LinearStabilityFactory<TModel>::createSolver()
{
   // Create model
   TModel model;
   model.init();

   // Create simulation
   auto spSolver = std::make_shared<MarginalCurve>();

   // Set spatial scheme
   auto spScheme = std::make_shared<typename TModel::SchemeType>(
      model.SchemeFormulation(), GridPurpose::SIMULATION);

   // Create list of field ID strings for boundary conditions
   std::vector<std::string> bcNames = model.backend().fieldNames();

   // Create list of nondimensional ID strings for physical parameters
   std::vector<std::string> ndNames = model.backend().paramNames();

   // Get model configuration tags
   auto modelCfg = model.configTags();

   // Add configuration file and parameters
   spSolver->setConfiguration(spScheme->dimension(), spScheme->tag(),
      model.backend().isPeriodicBox(), bcNames, ndNames, modelCfg);

   // Initialise simulation
   std::map<std::string, MHDFloat> cfg;
   std::set<SpatialScheme::Feature> features;
   spSolver->getConfig(cfg, features, model.version());
   spScheme->enable(features);
   model.configure(spScheme->features());
   spSolver->updateConfig(model.backend().automaticParameters(cfg));

   spSolver->initBase(model.spBackend());


   // Initialise resolution
   spSolver->initResolution(spScheme);

   // Set the boundary conditions
   SharedSimulationBoundary spBcs = spSolver->createBoundary();

   // Initialise the simulation
   spSolver->init(spBcs);

   return spSolver;
}

} // namespace QuICC

#endif // QUICC_LINEARSTABILITYFACTORY_HPP
