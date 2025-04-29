/**
 * @file RunStability.cpp
 * @brief General executable for the stability implementations
 */

/// Set the path to the simulation implementation
#define MODELPATH Model/QUICC_RUNSIM_PATH/PhysicalModel.hpp
/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
/// Define small macros allowing to convert to string
#define MAKE_STR( _P ) MAKE_STR_X( _P )
/// Create header include string for the required implementation
#define MODELHEADER MAKE_STR( MODELPATH )

// System includes
//

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Model/RunApplication.hpp"
#include "Stability/LinearStabilityFactory.hpp"
#include MODELHEADER

/**
 * @brief Setup and run the simulation
 */
int run()
{
   typedef QuICC::Model::QUICC_RUNSIM_CPPMODEL::PhysicalModel Model;
   int status = QuICC::run_application<QuICC::LinearStabilityFactory,Model>();

   return status;
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

   // Initialise everything that can't be done inside a class
   QuICC::StageTimer::allowIo(QuICC::QuICCEnv().allowsIO());
   QuICC::Profiler::Initialize();

   // Compute simulation
   QuICC::Profiler::RegionStart<1>("Walltime");
   auto code = run();
   QuICC::Profiler::RegionStop<1>("Walltime");

   // Finalise everything that can't be done inside a class
   QuICC::Profiler::Finalize();

   return code;
}
