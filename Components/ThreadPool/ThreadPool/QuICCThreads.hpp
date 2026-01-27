/**
 * @file QuICCThreads.hpp
 * @brief Static threadpool for QuICC
 */

#ifndef QUICC_QUICCTHREADS_HPP
#define QUICC_QUICCTHREADS_HPP

#ifdef QUICC_USE_THREADPOOL
// System includes
//

// Project includes
//
#include "ThreadPool/ThreadPool.hpp"

namespace QuICC {

ThreadPool::ThreadPool<ThreadPool::boost_t>& QuICCThreads();

} // namespace QuICC

// Dummy implementation
#else
   namespace QuICC {

   inline void QuICCThreads() {};

   } // namespace QuICC
#endif //QUICC_USE_THREADPOOL

#endif // QUICC_QUICCTHREADS_HPP
