/**
 * @file ChebyshevTests.cpp
 * @brief Catch2 tests driver for the Chebyshev DenseSM operators
 */

#define CATCH_CONFIG_RUNNER

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "TestSuite/DenseSM/Chebyshev/LinearMap/TestArgs.hpp"
#include "Profiler/Interface.hpp"

namespace test = QuICC::TestSuite::DenseSM::Chebyshev::LinearMap;

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

   std::string testType = "";

   // Build a new parser on top of Catch's
   using namespace Catch::clara;
   auto cli
      = session.cli()
      | Opt( test::args().ulp, "ulp" )        // Add max ulp option
         ["--ulp"]
         ("Maximum acceptable ulp")
      | Opt( test::args().params, "id" )      // Add test id option
         ["--id"]
         ("Test id")
      | Opt( testType, "test type" )          // Add test type
         ["--type"]
         ("Test type: ")
      | Opt( test::args().dumpData )          // Add keep output data option
         ["--dumpData"]
         ("Write output data to file?");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   // Set test from command line
   if(testType != "")
   {
      test::args().setType(testType);
   }

   if(test::args().params.size() > 0)
   {
      test::args().useDefault = false;
   }

   returnCode = session.run();

   QuICC::Profiler::Finalize();

   return returnCode;
}
