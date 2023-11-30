/**
 * @file GetVM.hpp
 * @brief Minimal utility to get process virtual memory usage.
 */

#ifndef QUICC_PROFILER_TRACKER_GETVM_HPP
#define QUICC_PROFILER_TRACKER_GETVM_HPP

#include <cstddef>


#if defined(__linux__) // Linux

#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>

#elif defined(__APPLE__) // OSX


#endif


namespace QuICC {
namespace Profiler {

/**
 * @brief 
 * @return virtual memory size in bytes
 *
 */
inline std::size_t GetVM()
{

#if defined(__linux__) // Linux

    std::ifstream stat_stream("/proc/self/stat", std::ios_base::in); //get info from proc
    
    std::string pid, comm, state, ppid, pgrp, session, tty_nr,
        tpgid, flags, minflt, cminflt, majflt, cmajflt,
        utime, stime, cutime, cstime, priority, nice,
        O, itrealvalue, starttime;
    std::size_t vsize;
    std::size_t rss;
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
        >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
        >> utime >> stime >> cutime >> cstime >> priority >> nice
        >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
    stat_stream.close();
    return vsize;

#elif defined(__APPLE__) // Darwin


#endif // Error or unknown OS

    return {};
}


} // namespace Profiler
} // namespace QuICC


#endif // QUICC_PROFILER_TRACKER_GETVM_HPP
