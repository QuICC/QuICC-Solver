/** 
 * @file ProfilerBase.cpp
 * @brief Source of the implementation of a profiling timer
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/ProfilerBase.hpp"

// Project includes
//

namespace QuICC {

namespace Debug {

namespace Profiler {

   const unsigned int ProfilerBase::mActiveLvl = LVL3;

   std::map<BreakPoint,MHDFloat> ProfilerBase::mTimings = std::map<BreakPoint,MHDFloat>();

   void ProfilerBase::init()
   {
   }

   bool ProfilerBase::isActive(BreakPoint point)
   {
      return ((static_cast<int>(point) % mActiveLvl) == 0);
   }

   MHDFloat ProfilerBase::time(BreakPoint point)
   {
      return mTimings.at(point);
   }

   void ProfilerBase::update(BreakPoint point, MHDFloat time)
   {
      if(mTimings.count(point) == 0)
      {
         mTimings.insert(std::make_pair(point, 0.0));
      }

      // Increment measured time
      mTimings.at(point) += time;
   }

   size_t ProfilerBase::size()
   {
      return mTimings.size();
   }

   void ProfilerBase::getTimings(Array& ts)
   {
      ts.resize(ProfilerBase::size());
      int i = 0;
      for(auto it = mTimings.begin(); it != mTimings.end(); ++it)
      {
         ts(i) = ProfilerBase::time(it->first);
         i++;
      }
   }

   std::map<BreakPoint,std::string>  ProfilerBase::namePoints()
   {
      std::map<BreakPoint,std::string>  names;
      for(auto it = mTimings.begin(); it != mTimings.end(); ++it)
      {
         names.insert(std::make_pair(it->first, pointName(it->first)));
      }

      return names;
   }

   void ProfilerBase::reset()
   {
      mTimings.clear();
   }

}
}
}
