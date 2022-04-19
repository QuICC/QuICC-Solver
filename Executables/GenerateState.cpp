/**
 * @file GenerateState.cpp
 * @brief Simple general executable to generate a state file for a model 
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
#include <stdexcept>

// Project includes
//
#include "QuICC/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/Model/StateGeneratorFactory.hpp"
#include MODELHEADER

// Check Model and Framework compatiblity
#include "QuICC/Version/ModelChecker.h"

/**
 * @brief Setup and run the simulation
 */
int run()
{
   int status = 0;

   // Create simulation
   QuICC::SharedStateGenerator   spGen;

   // Exception handling during the initialisation part
   try
   {
      // Create state generator
      spGen = QuICC::StateGeneratorFactory<QuICC::Model::QUICC_RUNSIM_CPPMODEL::PhysicalModel>::createGenerator();
   }

   // If exception is thrown, finalise (close files) and return
   catch(std::logic_error& e)
   {
      std::cerr << e.what() << std::endl;

      status = -1;
   }

   if(status == 0)
   {
      // Run the simulation
      spGen->run();
   }

   // Cleanup and close file handles
   spGen->finalize();

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
   // Initilise everything that can't be done inside a class
   QuICC::QuICCEnv().init();
   QuICC::StageTimer::allowIo(QuICC::QuICCEnv().allowsIO());

   // Storage for the return code of run
   int code;

   // Compute simulation
   code = run();

   // Finalise everything that can't be done inside a class
   QuICC::QuICCEnv().finalize();

   return code;
}
