/**
 * @file Affinity.cpp
 * @brief Source of the CPU affinity utils for thread pool
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
#include "ThreadPool/details/Affinity.hpp"

namespace QuICC {

namespace ThreadPool {

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

} // namespace details
} // namespace ThreadPool
} // namespace QuICC
