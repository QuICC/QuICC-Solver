/**
 * @file Tracker.hpp
 * @brief Simple profiler to track time and memory usage
 *
 * Track time and memory usage for a region.
 * Print summary with min/max/avg across mpi ranks.
 * Optionally write data collected in hdf5 format with highFive if
 * the macro QUICC_PROFILE_NATIVE_WRITER_HIGHFIVE is set.
 * The macro QUICC_PROFILE_NATIVE_SAMPLE defines the number of data
 * points collected.
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

#ifndef QUICC_PROFILE_NATIVE_SAMPLE
#define QUICC_PROFILE_NATIVE_SAMPLE 100
#endif

namespace QuICC {
namespace Profiler {

/**
 * @brief Simple profiler
 *
 * Track time and memory usage for a region.
 *
 */
class Tracker
{
public:
    static constexpr std::uint32_t collectSize = QUICC_PROFILE_NATIVE_SAMPLE;

    /**
     * tuple to store the tracked data
     * refer to enum #tracking for the meaning
     *
     */
    using tracking_t = std::tuple<
        std::array<double, collectSize>, double, double, std::uint32_t, SteadyTimer, MemoryTracker>;
    /**
     * list of data contained in #tracking_t
     *
     */
    enum tracking {time, memory, memoryDelta, count, timTracker, memTracker};
    /**
     * @brief initialize profiler
     * @param comm mpi communicator
     */
#ifdef QUICC_MPI
    static void init(MPI_Comm comm = MPI_COMM_WORLD);
#else
    static void init();
#endif
    /**
     * @brief initialize a region
     * @param name name of the region to profile
     *
     * Optionally initialize a region as a separate step from start.
     * It results is a lower overhead for the first invocation of start.
     */
    static void initRegion(const std::string& name);
    /**
     * @brief start profiling a region
     * @param name name of the region to profile
     *
     * If a region is not initialized, it is also initialized.
     */
    static void start(const std::string& name);
    /**
     * @brief stop profiling a region
     * @param name name of the region to profile
     */
    static void stop(const std::string& name);
    /**
     * @brief stop profiling a region
     * @param name name of the region to profile
     * @return tuple containing the data for a region
     */
    static tracking_t get(const std::string& name);
    /**
     * @brief print summary to a stream
     * @param os stream to use as output
     */
    static void print(std::ostream& os = std::cout);
    /**
     * @brief print all collected data to file with hdf5 format
     */
    static void print_hdf5();

private:

    /**
     * @brief map between region name and tracked data
     */
    static std::map<std::string, tracking_t> mRegions;

#ifdef QUICC_MPI
    /**
     * @brief copy of the mpi communicator
     */
    static MPI_Comm mComm;
#endif

};

} // namespace Profiler
} // namespace QuiCC
