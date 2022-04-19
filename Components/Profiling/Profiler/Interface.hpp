/**
 * @file Interface.hpp
 * @brief Source of the implementation of Profiler wrappers
 */

#include <string>
#if defined QUICC_MPI
#include <regex>
#include <mpi.h>
#endif

#if defined QUICC_PROFILER_BACKEND_NATIVE
#include "Tracker/Tracker.hpp"
#elif defined QUICC_PROFILER_BACKEND_LIKWID
#include <likwid.h>
#endif

namespace QuICC {

namespace Profiler {

namespace details {
    inline static void Initialize();
    inline static void Finalize();
    inline static void RegionStart(const std::string& name);
    inline static void RegionEnd(const std::string& name);

    class Barrier {
    public:
        static void Initialize()
        {
            #ifdef QUICC_MPI
            // Get enviroment variable to add mpi barriers
            const char* env = std::getenv("QUICC_PROFILER_MPI_BARRIER");
            if (env) {
                auto const regex_a = std::regex("after", std::regex::icase);
                barrier_after = std::regex_search(env, regex_a);
                auto const regex_b = std::regex("before", std::regex::icase);
                barrier_before = std::regex_search(env, regex_b);
            }
            #endif
        }
        inline static void Before()
        {
            #ifdef QUICC_MPI
            if (barrier_before) {
                MPI_Barrier(MPI_COMM_WORLD);
            }
            #endif
        }
        inline static void After()
        {
            #ifdef QUICC_MPI
            if (barrier_after) {
                MPI_Barrier(MPI_COMM_WORLD);
            }
            #endif
        }
    private:
        static bool barrier_before;
        static bool barrier_after;
    };

    // static members init
    bool Barrier::barrier_before = false;
    bool Barrier::barrier_after = false;

}


inline static void Initialize()
{
    details::Barrier::Initialize();
    details::Initialize();
}

inline static void Finalize()
{
    details::Finalize();
}

template <int L = 0>
inline static void RegionStart(const std::string& name)
{
    if constexpr (L <= QUICC_PROFILE_LEVEL)
    {
        details::Barrier::Before();
        details::RegionStart(name);
    }
}

template <int L = 0>
inline static void RegionEnd(const std::string& name)
{
    if constexpr (L <= QUICC_PROFILE_LEVEL)
    {
        details::RegionEnd(name);
        details::Barrier::After();
    }
}

// Provide RAII for Region
template <int L = 0>
class RegionFixture{
public:
    RegionFixture() = delete;
    RegionFixture(const std::string& InitName)
        : name(InitName) {
        RegionStart<L>(name);
    }
    ~RegionFixture() {
        RegionEnd<L>(name);
    }
private:
    const std::string& name;
};

//
//
// Specific implementations
//
//

#if defined QUICC_PROFILER_BACKEND_TESTER

namespace details {
    inline static void Initialize()
    {
        std::cout << "initall" << std::endl;
    }

    inline static void Finalize()
    {
        std::cout << "endall" << std::endl;
    }

    inline static void RegionStart(const std::string& name)
    {
        std::cout << "start: "+name << std::endl;
    }
    inline static void RegionEnd(const std::string& name)
    {
        std::cout << "end: "+name << std::endl;
    }
}

#elif defined QUICC_PROFILER_BACKEND_NATIVE

namespace details {
    inline static void Initialize()
    {
        Tracker::init();
    }

    inline static void Finalize()
    {
        Tracker::print();
    }

    inline static void RegionStart(const std::string& name)
    {
        Tracker::start(name);
    }
    inline static void RegionEnd(const std::string& name)
    {
        Tracker::stop(name);
    }
}

#elif defined QUICC_PROFILER_BACKEND_LIKWID

namespace details {
    inline static void Initialize()
    {
        likwid_markerInit();
    }

    inline static void Finalize()
    {
        likwid_markerClose();
    }

    inline static void RegionStart(const std::string& name)
    {
        likwid_markerStartRegion(name.c_str());
    }
    inline static void RegionEnd(const std::string& name)
    {
        likwid_markerStopRegion(name.c_str());
    }
}

#endif


} // namespace Profiler
} // namespace QuICC
