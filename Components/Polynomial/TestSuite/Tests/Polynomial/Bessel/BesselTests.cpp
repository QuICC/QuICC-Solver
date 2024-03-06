/**
 * @file BesselTests.cpp
 * @brief Catch2 tests driver for the Bessel operators
 */

#define CATCH_CONFIG_RUNNER

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Polynomial/Bessel/TestArgs.hpp"
namespace test = QuICC::TestSuite::Polynomial::Bessel;

int main(int argc, char* argv[])
{
   Catch::Session session; // There must be exactly one instance

   std::string testType = "";

   // Build a new parser on top of Catch's
   using namespace Catch::clara;
   auto cli =
      session.cli() |
      Opt(test::args().specN, "specN") // Add spectral truncation option
         ["--specN"]("Spectral truncation") |
      Opt(test::args().physN, "physN") // Add physical grid size option
         ["--physN"]("Physical grid size") |
      Opt(test::args().ulp, "ulp") // Add max ulp option
         ["--ulp"]("Maximum acceptable ulp") |
      Opt(test::args().params, "l") // Add harmonic degree option
         ["--l"]("Harmonic degree L") |
      Opt(test::args().ids, "test metadata ID") // Add test ID
         ["--id"]("Test metadata ID") |
      Opt(testType, "test type") // Add test type
         ["--type"](
            "Test type: matrix, weighted, otf_innner, otf_outer, otf_reduce") |
      Opt(test::args().dumpData) // Add dump output data option
         ["--dumpData"]("Write output data to file?");

   // Now pass the new composite back to Catch so it uses that
   session.cli(cli);

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine(argc, argv);
   if (returnCode != 0) // Indicates a command line error
   {
      return returnCode;
   }

   // Set test from command line
   if (testType != "")
   {
      test::args().setType(testType);
   }

   if (test::args().params.size() > 0)
   {
      if ((test::args().specN > 0 && test::args().physN > 0))
      {
         test::args().useDefault = false;
         test::args().ids.clear();
      }
      else
      {
         std::cerr << "You need to specify --l, --specN and --physN or provide "
                      "metadata ID with --id"
                   << std::endl;
         return 1;
      }
   }
   else if (test::args().ids.size() > 0)
   {
      test::args().useDefault = false;
   }

   return session.run();
}
