/**
 * @file SerialTimer.hpp
 * @brief Implementation of a serial timer 
 */

#ifndef QUICC_SERIALTIMER_HPP
#define QUICC_SERIALTIMER_HPP

// System includes
//
#include <ctime>

// External includes
//

// Project includes
//
#include "ITimer.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a serial timer
    */
   class SerialTimer: public ITimer
   {
      public:
         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         explicit SerialTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         virtual ~SerialTimer();

         /**
          * @brief Start clock
          */
         virtual void start();

         /**
          * @brief Stop clock
          */
         virtual void stop();

         /**
          * @brief Get elapsed time
          */
         virtual double time() const;

         /**
          * @brief Query current time without changes to state
          */
         virtual double queryTime() const;

         /**
          * @brief Reset timer (stop and restart)
          */
         virtual double resetTimer();
         
      protected:

      private:
         /**
          * @brief Start
          */
         timespec mStart;

         /**
          * @brief Stop
          */
         timespec mStop;

         /**
          * @brief Compute elapsed seconds between timer start and given stop
          */
         double elapsedSeconds(const timespec& stop) const;

         /**
          * @brief Compute elapsed seconds between timer start and timer stop
          */
         double elapsedSeconds() const;
   };

}

#endif // QUICC_SERIALTIMER_HPP
