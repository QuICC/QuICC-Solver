/**
 * @file ThreadPoolBoost.cpp
 * @brief Source of the Boost based thread pool
 */

// System includes
//
#include <thread>
#include <sched.h>
#include <boost/interprocess/detail/os_thread_functions.hpp>
#include <boost/thread/barrier.hpp>

// Project includes
//
#include "ThreadPool/ThreadPoolBoost.hpp"
#include "ThreadPool/details/Affinity.hpp"

namespace QuICC {

namespace ThreadPool {

ThreadPool<boost_t>::ThreadPool(const int size) : mSize(size)
{
   this->mspPool = std::make_unique<boost::asio::thread_pool>(size);

   // Configure affinity
   if(this->mSize > 0)
   {
      this->configureAffinity();
   }
}

int ThreadPool<boost_t>::size() const
{
   return this->mSize;
}

boost::asio::thread_pool& ThreadPool<boost_t>::pool()
{
   return *this->mspPool;
}

void ThreadPool<boost_t>::configureAffinity()
{
   // Don't overlap main thread with pool if enough cores are available
   auto ncores = details::countAffinity();
   int shift = 0;
   if(ncores > this->mSize)
   {
      shift = 1;
   }
   else if(ncores < this->mSize)
   {
      throw std::logic_error("Not enough cores available for thread pool");
   }

   // Set affinity of threads
   std::vector<std::future<int>> tasks;
   boost::barrier poolBar(this->mSize);
   for(int i = 0; i < this->mSize; i++)
   {
      int id = i + shift;
      auto fut = boost::asio::post(this->pool(), std::packaged_task<int()>(std::bind(details::setAffinity, std::ref(poolBar), id)));
      tasks.push_back(std::move(fut));
   }

   // Wait for threads
   bool failed = false;
   for(auto&& task: tasks)
   {
      task.wait();
      auto s = task.get();
      if(s == -1)
      {
         failed = true;
      }
   }

   if(failed)
   {
      throw std::logic_error("Failed to set affinity of thread pool");
   }

   // Set affinity for main thread
   boost::barrier mainBar(1);
   details::setAffinity(mainBar, 0);

//#ifdef QUICC_DEBUG
   details::printAffinity();
   tasks.clear();
   for(int i = 0; i < this->mSize; i++)
   {
      auto fut = boost::asio::post(this->pool(), std::packaged_task<int()>(&details::printAffinity));
      tasks.push_back(std::move(fut));
   }

   for(auto&& task: tasks)
   {
      task.wait();
   }
//#endif
}

} // namespace ThreadPool
} // namespace QuICC
