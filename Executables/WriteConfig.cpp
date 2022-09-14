/**
 * @file WriteConfig.cpp
 * @brief Simple executable to write a configuration file template for a model 
 */

/// Set the path to the simulation implementation
#define MODELPATH QuICC/Model/QUICC_RUNSIM_PATH/PhysicalModel.hpp
/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
/// Define small macros allowing to convert to string
#define MAKE_STR( _P ) MAKE_STR_X( _P )
/// Create header include string for the required implementation
#define MODELHEADER MAKE_STR( MODELPATH )

// Configuration includes
//

// System includes
//
#include <iostream>

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Io/Config/ConfigurationWriter.hpp"
#include "QuICC/Io/Config/Simulation/Physical.hpp"
#include "QuICC/Io/Config/Simulation/Boundary.hpp"
#include MODELHEADER

typedef QuICC::Model::QUICC_RUNSIM_CPPMODEL::PhysicalModel PModel;

/**
 * @brief Setup and run the simulation
 */
int run()
{
   // create model
   PModel model;
   model.init();

   // Create spatial scheme
   auto spScheme = std::make_shared<typename PModel::SchemeType>(model.SchemeFormulation(), QuICC::GridPurpose::SIMULATION);

   // Set dimension
   int dim = spScheme->dimension();

   // Set type string
   std::string type = spScheme->tag();

   // Box periodicity
   std::vector<bool> isPeriodicBox = model.backend().isPeriodicBox();

   // Create configuration writer
   QuICC::Io::Config::ConfigurationWriter writer(dim, isPeriodicBox, type);

   // Create list of field ID strings for boundary conditions
   std::vector<std::string> bcNames = model.backend().fieldNames();

   // Create list of nondimensional ID strings for physical parameters
   std::vector<std::string> ndNames = model.backend().paramNames();

   // Add the physical part
   auto spPhys = std::make_shared<QuICC::Io::Config::Simulation::Physical>(ndNames);
   writer.rspSimulation()->addNode(QuICC::Io::Config::Simulation::PHYSICAL, spPhys);

   // Add the boundary part
   auto spBound = std::make_shared<QuICC::Io::Config::Simulation::Boundary>(bcNames);
   writer.rspSimulation()->addNode(QuICC::Io::Config::Simulation::BOUNDARY, spBound);

   // Get model configuration tags
   auto modelCfg = model.configTags();

   // Add the model part
   writer.rspModel()->addNodes(modelCfg);

   // Initialise writer
   writer.init();

   // Write configuration
   writer.write();

   // Finalise writer
   writer.finalize();

   return 0;
}

/**
 * @brief General main, setting up MPI if required
 *
 * The actual program is in run to make sure MPI initialisations
 * are called before anything else and finalisation after destruction
 */
int main(int argc, char* argv[])
{
   // Environment fixture
   QuICC::QuICCEnv();

   // Initilise everything that can't be done inside a class
   QuICC::StageTimer::allowIo(QuICC::QuICCEnv().allowsIO());

   // Compute simulation
   auto code = run();

   return code;
}
