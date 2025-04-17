/**
 * @file QuICCThreads.cpp
 * @brief Source of the static threadpool for QuICC
 */

// Project includes
//
#include "ThreadPool/QuICCThreads.hpp"
#include "ThreadPool/ThreadPool.hpp"

namespace QuICC {

ThreadPool::ThreadPool<ThreadPool::boost_t>& QuICCThreads()
{
   static ThreadPool::ThreadPool<ThreadPool::boost_t> tp(ThreadPool::envSize(1));

   return tp;
}

} // namespace QuICC
