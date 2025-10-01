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
#include "QuICC/NonDimensional/registerStability.hpp"
#include "QuICC/NonDimensional/GrowthRate.hpp"
#include "QuICC/NonDimensional/MaxIteration.hpp"
#include "QuICC/NonDimensional/Nev.hpp"
#include "QuICC/NonDimensional/Omega.hpp"
#include "QuICC/NonDimensional/Sort.hpp"
#include "QuICC/NonDimensional/StabilityMode.hpp"
#include "QuICC/NonDimensional/Tolerance.hpp"
#include "QuICC/NonDimensional/WriteMtx.hpp"
#include "Stability/MarginalCurve.hpp"

namespace QuICC {

/**
 * @brief Implementation of the linear stability factory
 */
template <class TModel> class LinearStabilityFactory
{
public:
   /// Typedef of type of object created
   typedef std::shared_ptr<MarginalCurve> ReturnType;

   /**
    * @brief Create a shared simulation for the model
    */
   static ReturnType create();

   /**
    * @brief Constructor
    */
   LinearStabilityFactory() = delete;

   /**
    * @brief Destructor
    */
   ~LinearStabilityFactory() = delete;

protected:
   /**
    * @brief Add stability specific NonDimensional parameters
    */
   static void addParameters(std::vector<std::string>& ndNames);

private:
};

template <class TModel>
void LinearStabilityFactory<TModel>::addParameters(std::vector<std::string>& ndNames)
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
typename LinearStabilityFactory<TModel>::ReturnType LinearStabilityFactory<TModel>::create()
{
   // Register stability specific nondimensional parameters
   NonDimensional::registerStability();

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
   LinearStabilityFactory::addParameters(ndNames);

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
