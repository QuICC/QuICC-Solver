/**
 * @file RunApplication.hpp
 * @brief Generic driver for a QuICC application
 */

#ifndef QUICC_RUNAPPLICATION_HPP
#define QUICC_RUNAPPLICATION_HPP

// System includes
//
#include <iostream>
#include <stdexcept>

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#ifdef QUICC_USE_THREADPOOL
#include "ThreadPool/QuICCThreads.hpp"
#endif // QUICC_USE_THREADPOOL
#include "Profiler/Interface.hpp"
#include "QuICC/Timers/StageTimer.hpp"

namespace QuICC {

/**
 * @brief Setup and run application
 */
template <template <typename> class TFactory, typename TModel> int run_application()
{
   int status = 0;

   // Application pointer
   typename TFactory<TModel>::ReturnType   spApp;

   // Exception handling during the initialisation part
   try
   {
      // Create application
      spApp = TFactory<TModel>::create();
   }

   // If exception is thrown, finalise (close files) and return
   catch(std::logic_error& e)
   {
      try
      {
         QuICC::QuICCEnv().abort(e.what());
      }
      catch(std::logic_error& ee)
      {
         std::cerr << ee.what() << std::endl;
      }

      status = -1;
   }

   if(status == 0)
   {
#ifdef QUICC_USE_THREADPOOL
      // Initialize threads
      QuICC::QuICCThreads();
#endif // QUICC_USE_THREADPOOL

      // Run application
      spApp->run();

      // Cleanup and close file handles
      spApp->finalize();
   }

   return status;
}
}

#endif // QUICC_RUNAPPLICATION_HPP
