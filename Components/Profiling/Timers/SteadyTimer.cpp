/** 
 * @file SteadyTimer.cpp
 * @brief Source of the serial timer implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "SteadyTimer.hpp"

// Project includes
//

namespace QuICC {

   SteadyTimer::SteadyTimer(const bool autostart)
   {
      // Check if timer should be started at creation
      if(autostart)
      {
         this->start();
      }
   }

   SteadyTimer::~SteadyTimer()
   {
   }

   void SteadyTimer::start()
   {
      mInit = std::chrono::steady_clock::now();
   }

   void SteadyTimer::stop()
   {
      mFinal = std::chrono::steady_clock::now();
   }

   double SteadyTimer::time() const
   {
      // return elapsed seconds
      return this->elapsedSeconds();
   }

   double SteadyTimer::queryTime() const
   {
      // Get current stop
      auto stop = std::chrono::steady_clock::now();;

      double current = this->elapsedSeconds(stop);

      // return elapsed seconds
      return current;
   }

   double SteadyTimer::resetTimer()
   {
      // Stop the timer
      this->stop();

      // Get elapsed time
      double tmp = this->time();

      // Set start time to stopping time
      this->mInit = this->mFinal;

      // return elapsed time
      return tmp;
   }

} // namespace QuICC
