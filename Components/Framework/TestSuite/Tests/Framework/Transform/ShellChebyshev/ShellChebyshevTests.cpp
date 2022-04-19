/**
 * @file Initest.cpp
 * @brief Catch2 initialization of the Shell Chebyshev transform in Framework
 */

#define CATCH_CONFIG_RUNNER

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Framework/Transform/ShellChebyshev/TestArgs.hpp"

int main( int argc, char* argv[] )
{
   // Namespace alias for test
   namespace ns_test = QuICC::TestSuite::Framework::Transform::ShellChebyshev;
   // Test argument typedef
   typedef ns_test::TestArgs Args;

   Catch::Session session; // There must be exactly one instance

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
      | Opt( Args::blockSize, "blockSize" )    // Add transform block size option
         ["--blockSize"]
         ("Block size of transform")
      | Opt( Args::lower, "lower" )    // Add lower bound of physical grid option
         ["--lower"]
         ("Lower bound of physical grid")
      | Opt( Args::upper, "upper" )    // Add upper bound of physical grid option
         ["--upper"]
         ("Upper bound of physical grid")
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
