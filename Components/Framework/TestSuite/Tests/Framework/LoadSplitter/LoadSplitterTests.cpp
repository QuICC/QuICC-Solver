/**
 * @file LoadSplitterTests.cpp
 * @brief Catch2 tests driver for the LoadSplitter
 */

#define CATCH_CONFIG_RUNNER

// Configuration includes
//

// System includes
//
#include <catch2/catch.hpp>

// Project includes
//
#include "QuICC/TestSuite/Framework/LoadSplitter/TestArgs.hpp"
namespace test = QuICC::TestSuite::Framework::LoadSplitter;

int main( int argc, char* argv[] )
{
   Catch::Session session; // There must be exactly one instance

   // Build a new parser on top of Catch's
   using namespace Catch::clara;
   auto cli
      = session.cli()
      | Opt( test::args().op, "op" )      // Add test op option
         ["--op"]
         ("Test operator name")
      | Opt( test::args().algorithm, "algorithm" ) // Add test algorithm option
         ["--algorithm"]
         ("Splitting algorithm")
      | Opt( test::args().truncation, "" )      // Add test truncation scheme option
         ["--truncation"]
         ("Truncation scheme")
      | Opt( test::args().np, "np" )      // Add test np option
         ["--np"]
         ("# MPI ranks")
      | Opt( test::args().db, "db" )      // Add test db option
         ["--db"]
         ("ID of database file")
      | Opt( test::args().stage, "stage" )      // Add test stage option
         ["--stage"]
         ("ID of MPI stage")
      | Opt( test::args().dim1D, "dim1D" )      // Add test dim1D option
         ["--dim1d"]
         ("Test dim1D")
      | Opt( test::args().dim2D, "dim2D" )      // Add test dim2D option
         ["--dim2d"]
         ("Test dim2D")
      | Opt( test::args().dim3D, "dim3D" )      // Add test dim3D option
         ["--dim3d"]
         ("Test dim3D")
      | Opt( test::args().params, "id" )      // Add test id option
         ["--id"]
         ("MPI rank for which to generate splitting data")
      | Opt( test::args().factors, "factor" )   // Add imposed factor option
         ["--factor"]
         ("Factor to use in decomposition")
      | Opt( test::args().checkRanks )          // Add checkRanks option
         ["--checkRanks"]
         ("Check ranks individually?")
      | Opt( test::args().dumpData )          // Add dumpData option
         ["--dumpData"]
         ("Write output data to file?")
      | Opt( test::args().dumpDetails )       // Add dumpDetails option
         ["--dumpDetails"]
         ("Write detailed output data to file?");

   // Now pass the new composite back to Catch so it uses that
   session.cli( cli );

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   if(test::args().params.size() > 0)
   {
      test::args().useDefault = false;
   }

   return session.run();
}
