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
#include <iostream>

// Project includes
//
#include "ThreadPool/ThreadPoolBoost.hpp"

namespace QuICC {

namespace ThreadPool {

namespace details {
   int printAffinity();
   int countAffinity();
   int setAffinity(boost::barrier& bar, const int idx);
   }

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

#if 1
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
#endif
}

namespace details {

int printAffinity()
{
   cpu_set_t mask;
   int status = sched_getaffinity(0, sizeof(cpu_set_t), &mask);

   if(status == 0)
   {
      std::stringstream stream;
      stream << "Affinity of thread: " << std::this_thread::get_id() << std::endl;

      std::size_t nproc = boost::interprocess::ipcdetail::get_num_cores();

      stream << "sched_getaffinity = ";
      for(std::size_t i = 0; i < nproc; i++)
      {
         bool isSet = CPU_ISSET(i, &mask);
         stream << isSet << " ";
      }
      stream << std::endl;

      std::cerr << stream.str();
   }
   else
   {
      std::perror("sched_getaffinity");
   }

   return status;
}

int countAffinity()
{
   cpu_set_t mask;

   int status = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
   if(status == -1)
   {
      std::perror("sched_getaffinity");
      throw std::logic_error("Failed to get affinity mask");
   }
   std::size_t nproc = boost::interprocess::ipcdetail::get_num_cores();

   int count = 0;
   for(std::size_t i = 0; i < nproc; i++)
   {
      bool isSet = CPU_ISSET(i, &mask);
      if(isSet)
      {
         count++;
      }
   }

   return count;
}

int setAffinity(boost::barrier& bar, const int idx)
{
   cpu_set_t mask;

   // Get current affinity
   int status = sched_getaffinity(0, sizeof(cpu_set_t), &mask);
   if(status == 0)
   {
      std::size_t nproc = boost::interprocess::ipcdetail::get_num_cores();

      // Find CPU ID
      status = -1;
      int myCpu = -1;
      int cpuCount = -(idx + 1);
      for(std::size_t i = 0; i < nproc; i++)
      {
         bool isSet = CPU_ISSET(i, &mask);
         if(isSet)
         {
            cpuCount++;
            if(cpuCount == 0)
            {
               myCpu = i;
               status = 0;
               break;
            }
         }
      }
      if(status == 0)
      {
         CPU_ZERO(&mask);
         CPU_SET(myCpu, &mask);
         status = sched_setaffinity(0, sizeof(cpu_set_t), &mask);
         if(status == -1)
         {
            std::perror("sched_setaffinity");
         }
      }
   }
   else
   {
      std::perror("sched_getaffinity");
   }

   bar.wait();

   return status;
};

}

} // namespace ThreadPool
} // namespace QuICC
