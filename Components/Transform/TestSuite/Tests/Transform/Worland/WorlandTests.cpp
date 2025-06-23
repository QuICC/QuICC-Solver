/**
 * @file WorlandTests.cpp
 * @brief Catch2 tests driver for the Worland transforms
 */

#define CATCH_CONFIG_RUNNER

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/TestSuite/Transform/Worland/TestArgs.hpp"
#include "Profiler/Interface.hpp"

namespace test = QuICC::TestSuite::Transform::Worland;

int main( int argc, char* argv[] )
{
   QuICC::QuICCEnv();

   QuICC::Profiler::Initialize();

   Catch::Session session; // There must be exactly one instance

   std::string testType = "";

   std::string commandFile = "";

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
      | Opt( test::args().np, "np" )      // Add test np option
         ["--np"]
         ("# MPI ranks")
      | Opt( test::args().rank, "rank" )      // Add test rank option
         ["--rank"]
         ("MPI rank")
      | Opt( testType, "test type" )                          // Add test type
         ["--type"]
         ("Test type: projector, integrator, reductor, bfloop")
      | Opt( test::args().timeOnly )         // Add timing only
         ["--timeOnly"]
         ("Only time execution, don't check results")
      | Opt( test::args().iter, "iter" )     // Number of iterations
         ["--iter"]
         ("Iterations")
      | Opt( test::args().dumpData )          // Add keep output data option
         ["--dumpData"]
         ("Write output data to file?")
      | Opt( options_file )          // Read options from file
         ["--options_file"]
         ("Read command options from file");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   std::string testType = "";
   test::args().clear();

   std::vector<std::string> commands;
   QuICC::TestSuite::readLines(commands, options_file);

   for(const auto& command: commands)
   {
	   std::vector<std::string> options;
	   QuICC::TestSuite::splitLine(options, command, ';');

	   int run_argc;
	   std::vector<char *> run_argv;
   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( run_argc, run_argv.data() );
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

   auto ret = session.run();
   }

   QuICC::Profiler::Finalize();

   return ret;
}
