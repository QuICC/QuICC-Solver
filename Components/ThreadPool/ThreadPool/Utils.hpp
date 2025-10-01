/**
 * @file ThreadPool.hpp
 * @brief Utils for working with thread pool 
 */

#ifndef QUICC_THREADPOOL_UTILS_HPP
#define QUICC_THREADPOOL_UTILS_HPP

namespace QuICC {

namespace ThreadPool {

/**
 * @brief Get threadpool size from ENV variable
 */
int envSize(const int defaultSize = -1);

} // namespace ThreadPool
} // namespace QuICC


#endif // QUICC_THREADPOOL_UTILS_HPP
