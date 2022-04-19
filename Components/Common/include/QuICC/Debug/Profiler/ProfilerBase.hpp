/**
 * @file ProfilerBase.hpp
 * @brief Implementation of the base of a profiling timer 
 */

#ifndef QUICC_DEBUG_PROFILERBASE_HPP
#define QUICC_DEBUG_PROFILERBASE_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Debug/Profiler/BreakPoint.hpp"

namespace QuICC {

namespace Debug {

namespace Profiler {

   /**
    * @brief Implementation of the base of a profiling timer
    */
   class ProfilerBase
   {
      public:
         /**
          * @brief Initialise the profiler base
          */
         static void init();

         /**
          * @brief Get elapsed time for provided break point
          *
          * @param point Break point
          */
         static MHDFloat time(BreakPoint point);

         /**
          * @brief Reset the profiling timings
          */
         static void reset();

         /**
          * @brief Check if breakpoint is active
          *
          * @param point Break point
          */
         static bool isActive(BreakPoint point);

         /**
          * @brief Get map of used points to name
          */
         static std::map<BreakPoint,std::string> namePoints();

         /**
          * @brief Get number of breakpoints
          */
         static size_t size();

         /**
          * @brief Get the measured times
          *
          * @param ts   Timings for the breakpoints
          */
         static void getTimings(Array& ts);
         
      protected:
         /**
          * @brief Constructor
          */
         ProfilerBase();

         /**
          * @brief Destructor
          */
         virtual ~ProfilerBase();

         /**
          * @brief Update measured time
          *
          * @param point   Location that has been profiled
          * @param time    The measured time
          */
         static void update(BreakPoint point, MHDFloat time);

      private:
         /**
          * @brief Profiling timings
          */
         static const unsigned int mActiveLvl;

         /**
          * @brief Profiling timings
          */
         static std::map<BreakPoint,MHDFloat> mTimings;
   };

}
}
}

#endif // QUICC_DEBUG_PROFILERBASE_HPP
