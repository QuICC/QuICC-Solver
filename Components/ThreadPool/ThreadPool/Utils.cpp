/**
 * @file Utils.cpp
 * @brief Source of thread pool utils
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "ThreadPool/Utils.hpp"

namespace QuICC {

namespace ThreadPool {

int envSize(const int defaultSize)
{
   int size = defaultSize;

   const char* envSize = std::getenv("QUICC_THREADPOOL_SIZE");
   if (envSize)
   {
      size = std::stoi(envSize);
   }

   if(size < 0)
   {
      throw std::logic_error("Thread pool size not defined");
   }

   return size;
}

} // namespace ThreadPool
} // namespace QuICC
