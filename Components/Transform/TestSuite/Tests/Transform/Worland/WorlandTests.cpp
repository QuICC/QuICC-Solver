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
#include "TestSuite/Io.hpp"
#include "QuICC/TestSuite/Transform/Worland/TestArgs.hpp"
#include "Profiler/Interface.hpp"

namespace test = QuICC::TestSuite::Transform::Worland;

int main( int argc, char* argv[] )
{
   QuICC::QuICCEnv();

   QuICC::Profiler::Initialize();

   Catch::Session session; // There must be exactly one instance

   std::string testType = "";

   std::string optionsFile = "";

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
      | Opt( optionsFile, "options file" )          // Read options from file
         ["--options_file"]
         ("Read command options from file");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   if(optionsFile == "")
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
      std::vector<std::string> commands;
      QuICC::TestSuite::readLines(commands, optionsFile);

      std::vector<std::string> failed;
      std::size_t counter = 0;
      for(const auto& command: commands)
      {
         // Reset session
         QuICC::Profiler::RegionResetAll();
         Catch::ConfigData cd = {0};
         session.useConfigData(cd);
         testType = "";
         test::args().clear();

         // Process next line
         std::vector<std::string> options;
         QuICC::TestSuite::splitLine(options, command, ';');

         std::cout << "[" << counter << "/" << commands.size() << "]:";
         for(auto&& o: options)
         {
            std::cout << ' ' << o;
         }
         std::cout << std::endl;

         int run_argc;
         std::vector<char *> run_argv;
         QuICC::TestSuite::getCommand(run_argc, run_argv, options);

         // Let Catch (using Clara) parse the command line
         int ret = session.applyCommandLine( run_argc, run_argv.data() );
         if( ret != 0 ) // Indicates a command line error
            return ret;

         // Set test from command line
         if(testType != "")
         {
            test::args().setType(testType);
         }

         if(test::args().params.size() > 0)
         {
            test::args().useDefault = false;
         }

         ret = session.run();
         counter++;

         if(ret != 0)
         {
            failed.push_back(options.at(1));
         }
      }

      returnCode = failed.size();

      if(returnCode == 0)
      {
         std::cout << "All merged tests passed" << std::endl;
      }
      else
      {
         std::cout << "Following merged tests failed: " << std::endl;
         for(auto&& tag: failed)
         {
            std::cout << "   " << tag << std::endl;
         }
      }
   }

   QuICC::Profiler::Finalize();

   return returnCode;
}
