/**
 * @file QuICCThreads.hpp
 * @brief Static threadpool for QuICC
 */

#ifndef QUICC_QUICCTHREADS_HPP
#define QUICC_QUICCTHREADS_HPP

// System includes
//

// Project includes
//
#include "ThreadPool/ThreadPool.hpp"

namespace QuICC {

ThreadPool::ThreadPool<ThreadPool::boost_t>& QuICCThreads();

} // namespace QuICC

#endif // QUICC_QUICCTHREADS_HPP
