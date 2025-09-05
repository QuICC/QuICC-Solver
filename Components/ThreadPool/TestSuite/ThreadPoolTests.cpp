#include <catch2/catch.hpp>
#include <memory>
#include <iostream>
#include <boost/thread/thread.hpp>

// QuICC
#include "Profiler/Interface.hpp"
#include "ThreadPool/QuICCThreads.hpp"

TEST_CASE("Threadpool with void futures", "[void_futures]")
{
   auto& tp = QuICC::QuICCThreads();
   std::vector<std::future<void>> tasks;
   std::vector<std::size_t> data;
   std::mutex mtx;

   for(int i = 0; i < 10; i++)
   {
      data.emplace_back(0);
      auto fut = boost::asio::post(tp.pool(), std::packaged_task<void()>(
         [i, &data, &mtx]()
         {
            std::lock_guard<std::mutex> lock(mtx);
            std::cerr << "Task : " << i << " running in thread " << boost::this_thread::get_id() << std::endl;
         }
         )
      );
      tasks.push_back(std::move(fut));
   }

   // Wait for threads
   for(auto&& task: tasks)
   {
      task.wait();
   }

   // Check
   CHECK( true);
}

TEST_CASE("Threadpool with int futures", "[int_futures]")
{
   auto& tp = QuICC::QuICCThreads();
   std::vector<std::future<int>> tasks;

   // Check
   CHECK( true);
}
