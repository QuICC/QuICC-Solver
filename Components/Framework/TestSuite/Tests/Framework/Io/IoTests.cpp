/**
 * @file IoTests.cpp
 * @brief Catch2 tests driver for the IO tests
 */

#define CATCH_CONFIG_RUNNER

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/TestSuite/Framework/Io/TestArgs.hpp"
#include "Profiler/Interface.hpp"

namespace test = QuICC::TestSuite::Framework::Io;

int main( int argc, char* argv[] )
{
   // Environment fixture
   QuICC::QuICCEnv();

   #ifdef QUICC_MPI
   {
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      QuICC::QuICCEnv().setup(size);
   }
   #else
      QuICC::QuICCEnv().setup(1);
   #endif

   QuICC::Profiler::Initialize();

   Catch::Session session; // There must be exactly one instance

   // Build a new parser on top of Catch's
   using namespace Catch::clara;
   auto cli
      = session.cli()
      | Opt( test::args().ulp, "ulp" )       // Add max ulp option
         ["--ulp"]
         ("Maximum acceptable ulp")
      | Opt( test::args().params, "id" )     // Add test id option
         ["--id"]
         ("Test id")
      | Opt( test::args().timeOnly )         // Add timing only
         ["--timeOnly"]
         ("Only time execution, don't check results")
      | Opt( test::args().iter, "iter" )     // Number of iterations
         ["--iter"]
         ("Iterations")
      | Opt( test::args().dumpData )         // Add dumpData option
         ["--dumpData"]
         ("Write output data to file?");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   returnCode = session.run();

   QuICC::Profiler::Finalize();

   return returnCode;
}
