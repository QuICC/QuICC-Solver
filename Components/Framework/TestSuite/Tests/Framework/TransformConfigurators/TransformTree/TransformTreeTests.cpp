/**
 * @file TransformTreeTest.cpp
 * @brief Catch2 initialization of the transform tree in Framework
 */

#define CATCH_CONFIG_RUNNER

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformConfigurators/TransformTree/TestArgs.hpp"
namespace test = QuICC::TestSuite::Framework::TransformConfigurators::TransformTree;

int main( int argc, char* argv[] )
{
   typedef test::TestArgs Args;

   Catch::Session session; // There must be exactly one instance

   // Build a new parser on top of Catch's
   using namespace Catch::clara;
   auto cli
      = session.cli()
      | Opt( Args::keepData )          // Add keep output data option
         ["--keepData"]
         ("Write output data to file?");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   return session.run();
}
