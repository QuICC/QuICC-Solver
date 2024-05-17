
#define CATCH_CONFIG_RUNNER

#include <catch2/catch.hpp>
#include <memory>

#include "Profiler/Interface.hpp"

int main(int argc, char **argv)
{
  QuICC::Profiler::Initialize();

  Catch::Session session; // There must be exactly one instance

  auto returnCode = session.run();

  QuICC::Profiler::Finalize();

  return returnCode;
}
