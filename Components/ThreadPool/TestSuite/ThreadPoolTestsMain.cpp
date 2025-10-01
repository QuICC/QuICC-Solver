
#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <memory>

#include "Profiler/Interface.hpp"
#include "ThreadPool/QuICCThreads.hpp"

int main(int argc, char** argv)
{
   QuICC::QuICCThreads();
   QuICC::Profiler::Initialize();

   Catch::Session session; // There must be exactly one instance

   // Let Catch (using Clara) parse the command line
   int returnCode = session.applyCommandLine( argc, argv );
   if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

   returnCode = session.run();

   QuICC::Profiler::Finalize();

   return returnCode;
}
