/**
 * @file MpiProfiler.hpp
 * @brief Implementation of a MPI profiling timer 
 */

#ifndef QUICC_DEBUG_PROFILER_MPIPROFILER_HPP
#define QUICC_DEBUG_PROFILER_MPIPROFILER_HPP

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Debug/Profiler/ProfilerBase.hpp"

namespace QuICC {

namespace Debug {

namespace Profiler {

   /**
    * @brief Implementation of a MPI timer
    */
   class MpiProfiler: public ProfilerBase
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
          * @brief Stop clock
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
         MpiProfiler();

         /**
          * @brief Destructor
          */
         ~MpiProfiler();

      private:
         /**
          * @brief Start
          */
         static std::map<BreakPoint, MHDFloat> t_starts;

         /**
          * @brief Stop
          */
         static std::map<BreakPoint, MHDFloat> t_stops;
   };

}
}
}

#endif // QUICC_DEBUG_PROFILER_MPIPROFILER_HPP
