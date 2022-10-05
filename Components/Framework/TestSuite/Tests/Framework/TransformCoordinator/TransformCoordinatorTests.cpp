/**
 * @file TransformCoordinatorTests.cpp
 * @brief Catch2 tests driver for the TransformCoordinator
 */

#define CATCH_CONFIG_RUNNER

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Framework/TransformCoordinator/TestArgs.hpp"
namespace test = QuICC::TestSuite::Framework::TCoord;

int main( int argc, char* argv[] )
{
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
      | Opt( test::args().dim1D, "dim1D" )   // Add test dim1D option
         ["--dim1D"]
         ("Test dim1D")
      | Opt( test::args().dim2D, "dim2D" )   // Add test dim2D option
         ["--dim2D"]
         ("Test dim2D")
      | Opt( test::args().dim3D, "dim3D" )   // Add test dim3D option
         ["--dim3D"]
         ("Test dim3D")
      | Opt( test::args().algorithm, "" )    // Add test algorithm option
         ["--algorithm"]
         ("Comm splitting algorithm")
      | Opt( test::args().grouper, "" )      // Add test grouper option
         ["--grouper"]
         ("Comm grouping algorithm")
      | Opt( test::args().dumpData )         // Add dumpData option
         ["--dumpData"]
         ("Write output data to file?");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   if(test::args().dim1D > 0)
   {
      test::args().useDefault = false;
   }

   return session.run();
}
