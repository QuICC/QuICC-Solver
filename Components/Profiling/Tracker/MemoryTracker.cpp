/**
 * @file MemoryTracker.cpp
 * @brief Track maximum memory usage between 2 points
 */

#include "MemoryTracker.hpp"
#include "PeakRss.hpp"

namespace QuICC {
namespace Profiler {

void MemoryTracker::start(){
    mInit = PeakRss();
}

void MemoryTracker::stop(){
    mFinal = PeakRss();
}

/*
 * return value is in bytes
 */
std::size_t MemoryTracker::diff(){
    return mFinal - mInit;
}

} // namespace Profiler
} // namespace QuICC
