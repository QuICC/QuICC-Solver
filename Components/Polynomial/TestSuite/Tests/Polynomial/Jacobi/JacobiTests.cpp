/**
 * @file JacobiTests.cpp
 * @brief Catch2 tests driver for the Jacobi operators
 */

#define CATCH_CONFIG_RUNNER

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Polynomial/Jacobi/TestArgs.hpp"
namespace test = QuICC::TestSuite::Polynomial::Jacobi;

int main( int argc, char* argv[] )
{
   Catch::Session session; // There must be exactly one instance

   std::vector<double> alpha;
   std::vector<double> beta;
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
      | Opt( alpha, "jacobi_a" )              // Add jacobi alpha option
         ["--jacobi_a"]
         ("Jacobi alpha")
      | Opt( beta, "jacobi_b" )               // Add jacobi beta option
         ["--jacobi_b"]
         ("Jacobi beta")
      | Opt( testType, "test type" )                          // Add test type
         ["--type"]
         ("Test type: matrix, weighted, otf_innner, otf_outer, otf_reduce")
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

   if(alpha.size() != beta.size())
   {
      std::cerr << "Alpha and beta parameters need to be give in pairs" << std::endl;
      return 1;
   }
   else
   {
      for(int i = 0; i < alpha.size(); i++)
      {
         test::args().params.push_back(alpha.at(i));
         test::args().params.push_back(beta.at(i));
      }
   }

   if(test::args().params.size() > 0)
   {
      if(test::args().specN > 0 && test::args().physN > 0)
      {
         test::args().useDefault = false;
      } else
      {
         std::cerr << "You need to specify --specN and --physN" << std::endl;
         return 1;
      }
   }

   return session.run();
}
