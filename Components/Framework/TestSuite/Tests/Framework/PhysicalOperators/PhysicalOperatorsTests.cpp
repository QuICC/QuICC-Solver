 #define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>

#include "Environment/QuICCEnv.hpp"
#include "Profiler/Interface.hpp"

int main( int argc, char* argv[] )
{
   QuICC::QuICCEnv();
   QuICC::Profiler::Initialize();

   Catch::StringMaker<double>::precision = 15;

   int result = Catch::Session().run( argc, argv );

   QuICC::Profiler::Finalize();

   return result;
}
