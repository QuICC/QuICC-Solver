/**
 * @file SteadylTimer.hpp
 * @brief Implementation of a timer using steady::clock
 */

#ifndef QUICC_STEADYTIMER_HPP
#define QUICC_STEADYTIMER_HPP

// System includes
//
#include <chrono>

// External includes
//

// Project includes
//
#include "ITimer.hpp"

namespace QuICC {

   /**
    * @brief Implementation of a serial timer
    */
   class SteadyTimer: public ITimer
   {
      public:
         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         explicit SteadyTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         virtual ~SteadyTimer();

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

      private:
         /**
          * @brief Start
          */
         std::chrono::steady_clock::time_point mInit{};

         /**
          * @brief Stop
          */
         std::chrono::steady_clock::time_point mFinal{};

         /*
          * @brief return value is in seconds
          */
         inline double elapsedSeconds() const
         {
            return elapsedSeconds(mFinal);
         }
         inline double elapsedSeconds(const std::chrono::steady_clock::time_point& Final) const
         {
            return static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(
               Final - mInit).count()) * 1e-6;
         }
   };

} // namespace QuICC

#endif // QUICC_STEADYTIMER_HPP
