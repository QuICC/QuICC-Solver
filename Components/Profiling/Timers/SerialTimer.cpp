/** 
 * @file SerialTimer.cpp
 * @brief Source of the serial timer implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "SerialTimer.hpp"

// Project includes
//

namespace QuICC {

   SerialTimer::SerialTimer(const bool autostart)
   {
      // Check if timer should be started at creation
      if(autostart)
      {
         this->start();
      } else
      {
         // Zero initialise start timespec structures
         this->mStart.tv_sec = 0.0;
         this->mStart.tv_nsec = 0.0;

         // Zero initialise stop timespec structures
         this->mStop.tv_sec = 0.0;
         this->mStop.tv_nsec = 0.0;
      }
   }

   SerialTimer::~SerialTimer()
   {
   }

   void SerialTimer::start()
   {
      // Get the starting timespec
      clock_gettime(CLOCK_REALTIME, &this->mStart);
   }

   void SerialTimer::stop()
   {
      // Get the stopping timespec
      clock_gettime(CLOCK_REALTIME, &this->mStop);
   }

   double SerialTimer::time() const
   {
      // return elapsed seconds
      return this->elapsedSeconds();
   }

   double SerialTimer::queryTime() const
   {
      // Get current stop
      timespec stop;

      // Get current elasped time
      clock_gettime(CLOCK_REALTIME, &stop);
      double current = this->elapsedSeconds(stop);

      // return elapsed seconds
      return current;
   }

   double SerialTimer::resetTimer()
   {
      // Stop the timer
      this->stop();

      // Get elapsed time
      double tmp = this->time();

      // Set start time to stopping time
      this->mStart = this->mStop;

      // return elapsed time
      return tmp;
   }

   double SerialTimer::elapsedSeconds(const timespec& stop) const
   {
      // Compute elapsed seconds between the two timespecs
      return static_cast<double>(stop.tv_sec - this->mStart.tv_sec) + static_cast<double>(stop.tv_nsec - this->mStart.tv_nsec)/1.0e9;
   }

   double SerialTimer::elapsedSeconds() const
   {
      return this->elapsedSeconds(this->mStop);
   }

}
