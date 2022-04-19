/**
 * @file SerialProfiler.hpp
 * @brief Implementation of a serial profiling timer 
 */

#ifndef QUICC_DEBUG_PROFILER_SERIALPROFILER_HPP
#define QUICC_DEBUG_PROFILER_SERIALPROFILER_HPP

// System includes
//
#include <time.h>
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Debug/Profiler/ProfilerBase.hpp"

namespace QuICC {

namespace Debug {

namespace Profiler {

   /**
    * @brief Implementation of a serial timer
    */
   class SerialProfiler: public ProfilerBase
   {
      public:
         /**
          * @brief Initialise the timers
          */
         static void init();

         /**
          * @brief Reset the timers
          */
         static void reset();

         /**
          * @brief Start clock
          *
          * @param point Location that has been profiled
          */
         static void start(BreakPoint point);

         /**
          * @brief Stop clock and add time to storage
          *
          * @param point Location that has been profiled
          */
         static void stop(BreakPoint point);

         /**
          * @brief Analyze the measured times among whole framework
          *
          * @param ts   Timings for the breakpoints
          * @param min  Minimal value within framework
          * @param max  Maximal value within framework
          */
         static void analyze(Array& ts, Array& min, Array& max);
         
      protected:
         /**
          * @brief Constructor
          */
         SerialProfiler();

         /**
          * @brief Destructor
          */
         ~SerialProfiler();

      private:
         /**
          * @brief Starts
          */
         static std::map<BreakPoint, timespec> t_starts;

         /**
          * @brief Stops
          */
         static std::map<BreakPoint, timespec> t_stops;

         /**
          * @brief Compute elapsed seconds between start and stop
          */
         static MHDFloat elapsedSeconds(timespec &t1, timespec &t2);
   };

}
}
}

#endif // QUICC_DEBUG_PROFILER_SERIALPROFILER_HPP
