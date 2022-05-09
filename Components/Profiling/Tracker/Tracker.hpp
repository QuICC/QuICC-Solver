/*
 * Track time and memory usage for a region.
 * Print summary with min/max/avg across mpi ranks
 *
 */
#pragma once

#include "MemoryTracker.hpp"
#include "Timers/SteadyTimer.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#ifdef QUICC_MPI
#include <mpi.h>
#endif
#include <tuple>

namespace QuICC {
namespace Profiler {


class Tracker
{
public:

    using tracking_t = std::tuple<double, double, double, std::uint32_t, SteadyTimer, MemoryTracker>;
    enum tracking {time, memory, memoryDelta, count, timTracker, memTracker};

#ifdef QUICC_MPI
    static void init(MPI_Comm comm = MPI_COMM_WORLD);
#else
    static void init();
#endif
    static void initRegion(const std::string& name);
    static void start(const std::string& name);
    static void stop(const std::string& name);
    static tracking_t get(const std::string& name);
    static void print(std::ostream& os = std::cout);

private:

    // region name -> tracking tuple
    static std::map<std::string, tracking_t> mRegions;

#ifdef QUICC_MPI
    // MPI communicator
    static MPI_Comm mComm;
#endif

};

} // namespace Profiler
} // namespace QuiCC
