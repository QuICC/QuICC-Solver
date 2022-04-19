/**
 * @file Initest.cpp
 * @brief Catch2 initialization of the ALegendre transform in Framework
 */

#define CATCH_CONFIG_RUNNER

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Framework/Transform/ALegendre/TestArgs.hpp"
namespace test = QuICC::TestSuite::Framework::Transform::ALegendre;

int main( int argc, char* argv[] )
{
   typedef test::TestArgs Args;

   Catch::Session session; // There must be exactly one instance

   int harmonicM = -1;

   // Build a new parser on top of Catch's
   using namespace Catch::clara;
   auto cli
      = session.cli()
      | Opt( Args::specN, "specN" )    // Add spectral truncation option
         ["--specN"]
         ("Spectral truncation")
      | Opt( Args::physN, "physN" )    // Add physical grid size option
         ["--physN"]
         ("Physical grid size")
      | Opt( harmonicM, "harmonicM" )  // Add harmonic order option
         ["--harmonicM"]
         ("Harmonic order M")
      | Opt( Args::keepData )          // Add keep output data option
         ["--keepData"]
         ("Write output data to file?");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   if(harmonicM >= 0)
   {
      Args::harmonicM.push_back(harmonicM);

      if(Args::specN > 0 && Args::physN > 0)
      {
         Args::useDefault = false;
      }
   }

   return session.run();
}
