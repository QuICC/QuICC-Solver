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
         enum BreakPoint {INIT, PRERUN, FIRSTSTEP, RUN, POSTRUN, TOTAL, FIRST_RUN, FIRST_TIMESTEP, FIRST_NONLINEAR, TIMESTEP, NONLINEAR, NBREAKPOINT};

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
          * @brief Increment iteration
          */
         void iteration();

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
          * @brief Iteration counter
          */
         std::size_t mCount;

         /**
          * @brief Analyze the measured times
          */
         void analyze(Array& min, Array& max);

         /**
          * @brief Format timing for output
          */
         void formatTiming(std::ostream& stream, const std::string name, const BreakPoint pt, const Array& min, const Array& max, const int spaces, const int digits, const int base);
   };

}

#endif // QUICC_EXECUTIONTIMER_HPP
