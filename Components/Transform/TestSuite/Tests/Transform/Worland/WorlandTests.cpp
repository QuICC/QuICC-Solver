/**
 * @file WorlandTests.cpp
 * @brief Catch2 tests driver for the Worland transforms
 */

#define CATCH_CONFIG_RUNNER

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "TestSuite/Io.hpp"
#include "QuICC/TestSuite/Transform/Worland/TestArgs.hpp"
#include "QuICC/TestSuite/Transform/MergeTest.hpp"
#include "Profiler/Interface.hpp"

namespace test = QuICC::TestSuite::Transform::Worland;

int main( int argc, char* argv[] )
{
   QuICC::QuICCEnv();

   QuICC::Profiler::Initialize();

   Catch::Session session; // There must be exactly one instance

   std::string testType = "";

   QuICC::TestSuite::Transform::MergeInfo info;

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
      | Opt( test::args().np, "np" )         // Add test np option
         ["--np"]
         ("# MPI ranks")
      | Opt( test::args().rank, "rank" )     // Add test rank option
         ["--rank"]
         ("MPI rank")
      | Opt( testType, "test type" )         // Add test type
         ["--type"]
         ("Test type: projector, integrator, reductor, bfloop")
      | Opt( test::args().timeOnly )         // Add timing only
         ["--timeOnly"]
         ("Only time execution, don't check results")
      | Opt( test::args().iter, "iter" )     // Number of iterations
         ["--iter"]
         ("Iterations")
      | Opt( test::args().dumpData )         // Add keep output data option
         ["--dumpData"]
         ("Write output data to file?")
      | Opt( info.file, "options file" )     // Read options from file
         ["--options_file"]
         ("Read command options from file")
      | Opt( info.jid, "parallel id" )       // Read options from file
         ["--jid"]
         ("Id of parallel executor")
      | Opt( info.jN, "parallel jobs" )      // Read options from file
         ["--jN"]
         ("Number of parallel jobs");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
   {
      return returnCode;
   }

   // run with given command line arguments
   if(info.file == "")
   {
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
   }
   // Process commands from options file
   else
   {
      returnCode = QuICC::TestSuite::Transform::runMergedTests(info, session, test::args());
   }

   QuICC::Profiler::Finalize();

   return returnCode;
}
