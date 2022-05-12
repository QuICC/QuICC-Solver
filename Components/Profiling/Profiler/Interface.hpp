/**
 * @file Interface.hpp
 * @brief Source of the implementation of Profiler wrappers
 */
#ifndef QUICC_PROFILER_INTERFACE_HPP
#define QUICC_PROFILER_INTERFACE_HPP

#include <string>
#if defined QUICC_MPI
#include <regex>
#include <mpi.h>
#endif

#if defined QUICC_PROFILE_BACKEND_NATIVE
#include "Tracker/Tracker.hpp"
#elif defined QUICC_PROFILE_BACKEND_LIKWID
#include <likwid.h>
#endif

#if ! defined QUICC_PROFILE_LEVEL
#define QUICC_PROFILE_LEVEL -1
#endif

namespace QuICC {

namespace Profiler {

namespace details {
    inline static void Initialize();
    inline static void Finalize();
    inline static void RegionStart(const std::string& name);
    inline static void RegionStop(const std::string& name);

    class CheckMpi {
    public:
        static void Init()
        {
            #ifdef QUICC_MPI
            // Check that mpi is initialized
            int isMpiInit;
            MPI_Initialized(&isMpiInit);
            if(!isMpiInit)
            {
                MPI_Init(nullptr, nullptr);
                amItheOwner = true;
            }
            #endif
        }

        static void Finalize()
        {
            #ifdef QUICC_MPI
            if(amItheOwner)
            {
                MPI_Finalize();
            }
            #endif
        }
    private:
        inline static bool amItheOwner = false;
    };

    class Barrier {
    public:
        static void Initialize()
        {
            #ifdef QUICC_MPI
            // Get enviroment variable to add mpi barriers
            const char* env = std::getenv("QUICC_PROFILE_MPI_BARRIER");
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
        inline static bool barrier_before = false;
        inline static bool barrier_after = false;
    };
}


inline static void Initialize()
{
    details::CheckMpi::Init();
    details::Barrier::Initialize();
    details::Initialize();
}

inline static void Finalize()
{
    details::Finalize();
    details::CheckMpi::Finalize();
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
inline static void RegionStop(const std::string& name)
{
    if constexpr (L <= QUICC_PROFILE_LEVEL)
    {
        details::RegionStop(name);
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
        RegionStop<L>(name);
    }
private:
    const std::string name;
};

//
//
// Specific implementations
//
//

#if defined QUICC_PROFILE_BACKEND_TESTER

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
    inline static void RegionStop(const std::string& name)
    {
        std::cout << "stop: "+name << std::endl;
    }
}

#elif defined QUICC_PROFILE_BACKEND_NATIVE

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
    inline static void RegionStop(const std::string& name)
    {
        Tracker::stop(name);
    }
}

#elif defined QUICC_PROFILE_BACKEND_LIKWID

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
    inline static void RegionStop(const std::string& name)
    {
        likwid_markerStopRegion(name.c_str());
    }
}

#else // do nothing

namespace details {
    inline static void Initialize(){}
    inline static void Finalize(){}
    inline static void RegionStart(const std::string& name){}
    inline static void RegionStop(const std::string& name){}
}

#endif


} // namespace Profiler
} // namespace QuICC

#endif // QUICC_PROFILER_INTERFACE_HPP
