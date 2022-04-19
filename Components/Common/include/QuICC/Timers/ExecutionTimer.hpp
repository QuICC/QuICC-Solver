/**
 * @file ExecutionTimer.hpp
 * @brief Implementation of an execution timer 
 */

#ifndef QUICC_EXECUTIONTIMER_HPP
#define QUICC_EXECUTIONTIMER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "Timers/TimerMacro.h"

namespace QuICC {

   /**
    * @brief Implementation of an execution timer
    */
   class ExecutionTimer: public TimerMacro
   {
      public:
         /**
          * @brief Create list of possible timing intervals
          */
         enum BreakPoint {INIT, PRERUN, RUN, POSTRUN, TOTAL};

         /**
          * @brief Constructor
          *
          * @param autostart Should the timer start at creation ?
          */
         explicit ExecutionTimer(const bool autostart = false);

         /**
          * @brief Destructor
          */
         ~ExecutionTimer();

         /**
          * @brief Update the timing of given region. 
          *
          * The possible breakpoints are defined in the BreakPoint enum
          *
          * @param point Breakpoint for which the timing has to be updated
          */
         void update(BreakPoint point);

         /**
          * @brief Query current time of given region. 
          *
          * The possible breakpoints are defined in the BreakPoint enum
          *
          * @param point Breakpoint for which the timing has to be updated
          */
         double queryTime(BreakPoint point) const;
         
         /**
          * @brief Print execution time information to stream
          *
          * @param stream  Output stream
          */
         void printInfo(std::ostream& stream);

      protected:

      private:

         using TimerMacro::queryTime;

         /**
          * @brief Storage for the execution times
          */
         Array mTimes;

         /**
          * @brief Analyze the measured times
          */
         void analyze(Array& min, Array& max);
   };

}

#endif // QUICC_EXECUTIONTIMER_HPP
