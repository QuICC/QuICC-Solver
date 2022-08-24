/**
 * @file MemoryTracker.hpp
 * @brief Track maximum memory usage between 2 points
 */

#ifndef QUICC_PROFILER_TRACKER_MEMORYTRACKER_HPP
#define QUICC_PROFILER_TRACKER_MEMORYTRACKER_HPP

#include <cstddef>

namespace QuICC {
namespace Profiler {

/**
 * @brief Track maximum memory usage between 2 points
 *
 * Tracks memory usage (peak working set size) between 2 points.
 * Modeled after a classic timer.
 */
class MemoryTracker {

public:
    /**
     * @brief start tracking memory
     */
    void start();
    /**
     * @brief stop tracking memory
     */
    void stop();
    /**
     * @brief compute difference between start and stop
     * @return value is in bytes
     */
    std::size_t diff();
private:
    /**
     * @brief store value at start
     */
    std::size_t mInit{};
    /**
     * @brief store value at stop
     */
    std::size_t mFinal{};
};

} // namespace Profiler
} // namespace QuICC

#endif // QUICC_PROFILER_TRACKER_MEMORYTRACKER_HPP
