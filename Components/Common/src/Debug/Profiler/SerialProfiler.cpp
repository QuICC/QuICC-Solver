/** 
 * @file SerialProfiler.cpp
 * @brief Source of the serial profiler implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/SerialProfiler.hpp"

// Project includes
//

namespace QuICC {

namespace Debug {

namespace Profiler {

   std::map<BreakPoint, timespec> SerialProfiler::t_starts = std::map<BreakPoint, timespec>();

   std::map<BreakPoint, timespec> SerialProfiler::t_stops = std::map<BreakPoint, timespec>();

   void SerialProfiler::init()
   {
      ProfilerBase::init();
   }

   void SerialProfiler::reset()
   {
      t_starts.clear();
      t_stops.clear();

      ProfilerBase::reset();
   }

   void SerialProfiler::start(BreakPoint point)
   {
      if(isActive(point))
      {
         if(t_starts.count(point) == 0)
         {
            t_starts.insert(std::make_pair(point, timespec()));
            t_stops.insert(std::make_pair(point, timespec()));
         }

         // Get the starting timespec
         clock_gettime(CLOCK_REALTIME, &t_starts.at(point));
      }
   }

   void SerialProfiler::stop(BreakPoint point)
   {
      if(isActive(point))
      {
         // Get the stopping timespec
         clock_gettime(CLOCK_REALTIME, &t_stops.at(point));

         // Store time
         ProfilerBase::update(point, SerialProfiler::elapsedSeconds(t_starts.at(point), t_stops.at(point)));
      }
   }

   MHDFloat SerialProfiler::elapsedSeconds(timespec &t1, timespec &t2)
   {
      // Compute elapsed seconds between the two timespecs
      return static_cast<MHDFloat>(t2.tv_sec - t1.tv_sec) + static_cast<MHDFloat>(t2.tv_nsec - t1.tv_nsec)/1.0e9;
   }

   void SerialProfiler::analyze(Array& ts, Array& min, Array& max)
   {
      SerialProfiler::getTimings(ts);

      min.resize(0);
      max.resize(0);
   }

}
}
}
