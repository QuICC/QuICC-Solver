#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

#include <chrono>
#include <thread>

#include "Tracker/GetVM.hpp"

// save
#if defined QUICC_PROFILE_LEVEL
    #define GLOBAL_QUICC_PROFILE_LEVEL QUICC_PROFILE_LEVEL
    #undef QUICC_PROFILE_LEVEL
#endif
#define QUICC_PROFILE_LEVEL 1
#if ! defined QUICC_PROFILE_BACKEND_NATIVE
    #define QUICC_PROFILE_BACKEND_NATIVE
    #define QUICC_PROFILE_BACKEND_NATIVE_CLEANDEF
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

TEST_CASE("Open Close", "[OpenClose]")
{
    using namespace QuICC::Profiler;
    Initialize();
    Finalize();
}

TEST_CASE("Start Stop region with String", "[StartStopString]")
{
    using namespace QuICC::Profiler;
    using namespace std::chrono_literals;
    Initialize();
    std::string region = "Count times hundred millisec";
    std::size_t count = 10;
    for (std::size_t i = 0; i < count; ++i)
    {
        RegionStart(region);
            std::this_thread::sleep_for(100ms);
        RegionStop(region);
    }

    REQUIRE(std::get<Tracker::tracking::count>(Tracker::get(region)) == count);
    // This check cannot be very precise...
    // Sum up time
    double timeTot{};
    for (std::size_t i = 0; i < count; ++i)
    {
        timeTot += std::get<Tracker::tracking::time>(Tracker::get(region))[i];
    }
    double timeRef = 0.1*count;
    REQUIRE(std::abs(timeTot - timeRef) < timeRef/10);

    Finalize();
}


TEST_CASE("Nested Levels: on, on, off", "[Levels]")
{
    using namespace QuICC::Profiler;
    using namespace std::chrono_literals;
    std::size_t maxVM = 0;
    Initialize();
    {
        auto vm0 = GetVM();
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
        auto vm1 = GetVM();
        maxVM =  std::max(maxVM, vm1 - vm0);
    }
    REQUIRE(std::get<Tracker::tracking::count>(Tracker::get("Main")) == 1);
    REQUIRE(std::get<Tracker::tracking::count>(Tracker::get("Nested")) == 1);
    REQUIRE(std::get<Tracker::tracking::count>(Tracker::get("Innermost")) == 0);

    // The following checks cannot be very precise...
    REQUIRE(std::abs(std::get<Tracker::tracking::time>(Tracker::get("Main"))[0] - 0.3) < 0.1);
    REQUIRE(std::abs(std::get<Tracker::tracking::time>(Tracker::get("Nested"))[0] - 0.2) < 0.1);
    //
    auto maxVMMB = maxVM * 1.0e-6;
    REQUIRE(std::abs(maxVMMB - 80.0) < 0.02);
    REQUIRE(std::get<Tracker::tracking::memoryDelta>(Tracker::get("Main")) <= maxVMMB * 1.12);
    REQUIRE(std::abs(std::get<Tracker::tracking::memoryDelta>(Tracker::get("Nested")) - 0.0) < 0.01);

    // high watermark is hard to estimate a priori
    Finalize();
}


// restore
#undef QUICC_PROFILE_LEVEL
#if defined GLOBAL_QUICC_PROFILE_LEVEL
    #define QUICC_PROFILE_LEVEL GLOBAL_QUICC_PROFILE_LEVEL
    #undef GLOBAL_QUICC_PROFILE_LEVEL
#endif
#if defined QUICC_PROFILE_BACKEND_NATIVE_CLEANDEF
    #undef QUICC_PROFILE_BACKEND_NATIVE_CLEANDEF
    #undef QUICC_PROFILE_BACKEND_NATIVE
#endif
