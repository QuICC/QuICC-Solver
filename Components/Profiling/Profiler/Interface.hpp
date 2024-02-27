/**
 * @file Interface.hpp
 * @brief Implementation of profiler wrappers
 *
 * Implementation of profiler wrappers. A backend is selected by a macro
 * QUICC_PROFILE_BACKEND_<name>.
 * A new backend needs to implement 4 methods (Initialize, Finalize, RegionStart, RegionStop, RegionResetAll)
 * in the namespace QuICC::Profiler::details.
 * A region is profiler only if the macro QUICC_PROFILE_LEVEL is equal
 * or higher than the template parameter L that is used for the region.
 *
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
    /**
     * @brief Internal interface for backend initialization
     *
     * Every backend needs to implement this method.
     */
    inline static void Initialize();
    /**
     * @brief Internal interface for backend finalization
     *
     * Every backend needs to implement this method.
     */
    inline static void Finalize();
    /**
     * @brief Internal interface to start the profiler for a region
     * @param name name of the region to profile
     *
     * Every backend needs to implement this method.
     */
    inline static void RegionStart(const std::string& name);
    /**
     * @brief Internal interface to stop the profiler for a region
     * @param name name of the region to profile
     *
     * Every backend needs to implement this method.
     */
    inline static void RegionStop(const std::string& name);
    /**
     * @brief Internal interface to reset the profiler for all regions
     *
     * Every backend needs to implement this method.
     */
    inline static void RegionResetAll();

    /**
     * @brief Class to check wether mpi is already initialized
     *
     * This is needed internally to make sure mpi is initialized.
     * If not it is initialized and finalized internally.
     */
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

    /**
     * @brief Class to add barriers to profiled regions
     *
     * This is needed internally to add a barrier before, after (or both) a region.
     * The barrier is activated if an enviroment variable (QUICC_PROFILE_MPI_BARRIER) contains
     * the string <after> and/or <before.
     */
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

/**
 * @brief Public interface to initialize the profiler
 */
inline static void Initialize()
{
    details::CheckMpi::Init();
    details::Barrier::Initialize();
    details::Initialize();
}

/**
 * @brief Public interface to finalize the profiler
 */
inline static void Finalize()
{
    details::Finalize();
    details::CheckMpi::Finalize();
}

/**
 * @brief Public interface to start the profiler for a region
 * @param name name of the region to profile
 * @tparam L profiling level
 */
template <int L = 0>
inline static void RegionStart(const std::string& name)
{
    if constexpr (L <= QUICC_PROFILE_LEVEL)
    {
        details::Barrier::Before();
        details::RegionStart(name);
    }
}

/**
 * @brief Public interface to stop the profiler for a region
 * @param name name of the region to profile
 * @tparam L profiling level
 */
template <int L = 0>
inline static void RegionStop(const std::string& name)
{
    if constexpr (L <= QUICC_PROFILE_LEVEL)
    {
        details::RegionStop(name);
        details::Barrier::After();
    }
}

/**
 * @brief Public interface to reset the profiler counter for all regions
 */
inline static void RegionResetAll()
{
    details::RegionResetAll();
}

/**
 * @brief Provide RAII functionality for a region to be profiled
 * @tparam L profiling level
 */
template <int L = 0>
class RegionFixture{
public:
    RegionFixture() = delete;
    /**
     * @brief fixture ctor
     * @param InitName name of the region to profile
     */
    RegionFixture(const std::string& InitName)
        : name(InitName) {
        RegionStart<L>(name);
    }
    ~RegionFixture() {
        RegionStop<L>(name);
    }
private:
    /**
     *@brief storage for the name of the region to profile
     */
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
    inline static void RegionResetAll()
    {
        std::cout << "reset all regions " << std::endl;
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
        Tracker::print_hdf5();
    }

    inline static void RegionStart(const std::string& name)
    {
        Tracker::start(name);
    }
    inline static void RegionStop(const std::string& name)
    {
        Tracker::stop(name);
    }
    inline static void RegionResetAll()
    {
        Tracker::resetAll();
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
    inline static void RegionResetAll()
    {            
        // we need to call likwid_markerResetRegion(reg->second) on all the regions - how to get the list of them? 
    }
}

#else // do nothing

namespace details {
    inline static void Initialize(){}
    inline static void Finalize(){}
    inline static void RegionStart(const std::string& name){}
    inline static void RegionStop(const std::string& name){}
    inline static void RegionResetAll(){}
}

#endif


} // namespace Profiler
} // namespace QuICC

#endif // QUICC_PROFILER_INTERFACE_HPP
