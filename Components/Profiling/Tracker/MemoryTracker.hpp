/**
 * @file MemoryTracker.hpp
 * @brief Track maximum memory usage between 2 points
 */

#pragma once

#include <cstddef>

namespace QuICC {
namespace Profiler {


class MemoryTracker {

public:
    void start();
    void stop();

    /*
     * return value is in bytes
     */
    std::size_t diff();
private:
    std::size_t mInit{};
    std::size_t mFinal{};
};

} // namespace Profiler
} // namespace QuICC
