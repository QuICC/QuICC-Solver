/**
 * @file QuadratureTests.cpp
 * @brief Catch2 tests driver for the quadrature rules
 */

#define CATCH_CONFIG_RUNNER

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Polynomial/Quadrature/TestArgs.hpp"
namespace test = QuICC::TestSuite::Polynomial::Quadrature;

int main( int argc, char* argv[] )
{
   Catch::Session session; // There must be exactly one instance

   std::string testType = "";

   // Build a new parser on top of Catch's
   using namespace Catch::clara;
   auto cli
      = session.cli()
      | Opt( test::args().specN, "specN" )    // Add spectral truncation option
         ["--specN"]
         ("Spectral truncation")
      | Opt( test::args().physN, "physN" )    // Add physical grid size option
         ["--physN"]
         ("Physical grid size")
      | Opt( test::args().ulp, "ulp" )        // Add max ulp option
         ["--ulp"]
         ("Maximum acceptable ulp")
      | Opt( test::args().params, "harmonicL" )  // Add harmonic degree option
         ["--harmonicL"]
         ("Harmonic degree L")
      | Opt( test::args().dumpData )          // Add dump output data option
         ["--dumpData"]
         ("Write output data to file?");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   // Set quadrature test
   test::args().setType("quadrature");

   if(test::args().physN > 0)
   {
      test::args().useDefault = false;
      test::args().gridN.clear();
      test::args().gridN.push_back(test::args().physN);
   }

   return session.run();
}
