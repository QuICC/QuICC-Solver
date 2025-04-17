/**
 * @file ThreadPool.hpp
 * @brief Dispatch header for threadpool implementation
 */

#ifndef QUICC_THREADPOOL_THREADPOOL_HPP
#define QUICC_THREADPOOL_THREADPOOL_HPP

// Boost implementation of thread pool
#include "ThreadPool/ThreadPoolBoost.hpp"

namespace QuICC {

namespace ThreadPool {

/**
 * @brief Get threadpool size from ENV variable
 */
int envSize(const int defaultSize = -1);

} // namespace ThreadPool
} // namespace QuICC


#endif // QUICC_THREADPOOL_THREADPOOL_HPP
