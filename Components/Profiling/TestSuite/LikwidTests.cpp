#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

#include <chrono>
#include <thread>

// save
#if defined QUICC_PROFILE_LEVEL
    #define GLOBAL_QUICC_PROFILE_LEVEL QUICC_PROFILE_LEVEL
    #undef QUICC_PROFILE_LEVEL
#endif
#define QUICC_PROFILE_LEVEL 1
#if ! defined QUICC_PROFILER_BACKEND_LIKWID
    #define QUICC_PROFILER_BACKEND_LIKWID
    #define QUICC_PROFILER_BACKEND_LIKWID_CLEANDEF
#endif

#include "Profiler/Interface.hpp"


int main( int argc, char* argv[] ) {
    #ifdef QUICC_MPI
        MPI_Init(nullptr, nullptr);
    #endif

    int result = Catch::Session().run( argc, argv );

    #ifdef QUICC_MPI
        MPI_Finalize();
    #endif

    return result;
}

// This is only an integration test
TEST_CASE("Nested Levels: on, on, off", "[Levels]")
{
    using namespace QuICC::Profiler;
    using namespace std::chrono_literals;
    Initialize();
    {
        RegionFixture<0> mainFix("Main");
        std::this_thread::sleep_for(100ms);
        std::vector<double> page(1e7);
        {
            RegionFixture<1> nestedFix("Nested");
            std::this_thread::sleep_for(100ms);
            {
                RegionFixture<2> innerFix("Innermost");
                std::this_thread::sleep_for(100ms);
            }
        }
    }

    Finalize();
}

// restore
#undef QUICC_PROFILE_LEVEL
#if defined GLOBAL_QUICC_PROFILE_LEVEL
    #define QUICC_PROFILE_LEVEL GLOBAL_QUICC_PROFILE_LEVEL
    #undef GLOBAL_QUICC_PROFILE_LEVEL
#endif
#if defined QUICC_PROFILER_BACKEND_LIKWID_CLEANDEF
    #undef QUICC_PROFILER_BACKEND_LIKWID_CLEANDEF
    #undef QUICC_PROFILER_BACKEND_LIKWID
#endif
