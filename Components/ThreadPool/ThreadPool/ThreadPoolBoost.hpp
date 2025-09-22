/**
 * @file ThreadPoolBoost.hpp
 * @brief Thread pool using Boost
 */

#ifndef QUICC_THREADPOOL_THREADPOOLBOOST_HPP
#define QUICC_THREADPOOL_THREADPOOLBOOST_HPP

// System includes
//
#include <boost/asio.hpp>
#include <future>

// Project includes
//
#include "ThreadPool/Tags.hpp"

namespace QuICC {

namespace ThreadPool {

template <typename T> class ThreadPool;

/**
 * @brief Thread pool using Boost
 */
template <> class ThreadPool<boost_t>
{
public:
   /**
    * @brief Constructor
    */
   explicit ThreadPool(const int size);

   /**
    * @brief Constructor
    */
   virtual ~ThreadPool() = default;

   /**
    * @Brief Pool size
    */
   int size() const;

   /**
    * @brief Get thread pool
    */
   boost::asio::thread_pool& pool();

protected:
private:
   /**
    * @brief Configure thread's CPU affinity
    */
   void configureAffinity();

   /**
    * @brief Thread pool
    */
   std::unique_ptr<boost::asio::thread_pool> mspPool;

   /**
    * @brief Size
    */
   const int mSize;
};

} // namespace ThreadPool
} // namespace QuICC

#endif // QUICC_THREADPOOL_THREADPOOLBOOST_HPP
