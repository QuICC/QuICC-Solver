#include <chrono>
#include <thread>

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

int main()
{
    using namespace std::chrono_literals;
    MPI_Init(nullptr, nullptr);
    int rank{};
    QuICC::Profiler::Initialize();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    {
        QuICC::Profiler::RegionFixture mainFix("main");
        for (int i = 0; i < 10; ++i)
        {
            QuICC::Profiler::RegionFixture loopFix("loop");
            if (rank == 0)
            {
                std::this_thread::sleep_for(100ms);
            }
            else
            {
                std::vector<double> v(1e6, 1.0);
                std::this_thread::sleep_for(200ms);
                // avoid opt out
                v[2] = 3.0;
                std::cout << v[3];
            }
        }
        std::cout << std::endl;
    }
    QuICC::Profiler::Finalize();
    MPI_Finalize();

    return 0;
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
