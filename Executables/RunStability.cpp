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

// Configuration includes
//

// System includes
//
#include <iostream>
#include <stdexcept>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "Stability/LinearStabilityFactory.hpp"
#include "Stability/MarginalCurve.hpp"
#include MODELHEADER
#include "Profiler/Interface.hpp"

/**
 * @brief Setup and run the simulation
 */
int run()
{
   QuICC::Profiler::RegionFixture runFix("run");

   int status = 0;

   // Create simulation
   std::shared_ptr<QuICC::MarginalCurve> spSolver;

   // Exception handling during the initialisation part
   try
   {
      // Create simulation
      spSolver = QuICC::LinearStabilityFactory<
         QuICC::Model::QUICC_RUNSIM_CPPMODEL::PhysicalModel>::createSolver();
   }

   // If exception is thrown, finalise (close files) and return
   catch (std::logic_error& e)
   {
      try
      {
         QuICC::QuICCEnv().abort(e.what());
      }
      catch (std::logic_error& ee)
      {
         std::cerr << ee.what() << std::endl;
      }

      status = -1;
   }

   if (status == 0)
   {
      // Run the simulation
      spSolver->run();

      // Cleanup and close file handles
      spSolver->finalize();
   }

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
   auto code = run();

   // Finalise everything that can't be done inside a class
   QuICC::Profiler::Finalize();

   return code;
}
