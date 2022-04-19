/** 
 * @file MpiProfiler.cpp
 * @brief Source of the implementation of MPI profiler
 */

// Configuration includes
//

// System includes
//
#include <mpi.h>

// External includes
//

// Class include
//
#include "QuICC/Debug/Profiler/MpiProfiler.hpp"

// Project includes
//

namespace QuICC {

namespace Debug {

namespace Profiler {

   std::map<BreakPoint, MHDFloat> MpiProfiler::t_starts = std::map<BreakPoint, MHDFloat>();

   std::map<BreakPoint, MHDFloat> MpiProfiler::t_stops = std::map<BreakPoint, MHDFloat>();

   void MpiProfiler::init()
   {
      ProfilerBase::init();
   }

   void MpiProfiler::reset()
   {
      t_starts.clear();
      t_stops.clear();

      ProfilerBase::reset();
   }

   void MpiProfiler::start(BreakPoint point)
   {
      if(isActive(point))
      {
         if(t_starts.count(point) == 0)
         {
            t_starts.insert(std::make_pair(point, 0));
            t_stops.insert(std::make_pair(point, 0));
         }

         // Get starting MPI time
         t_starts.at(point) = MPI_Wtime();
      }
   }

   void MpiProfiler::stop(BreakPoint point)
   {
      if(isActive(point))
      {
         // Get stoppping MPI time
         t_stops.at(point) = MPI_Wtime();

         // Store time
         ProfilerBase::update(point, t_stops.at(point)-t_starts.at(point));
      }
   }

   void MpiProfiler::analyze(Array& ts, Array& min, Array& max)
   {
      MpiProfiler::getTimings(ts);

      // Resize the storage
      min.resize(ts.size());
      max.resize(ts.size());

      // Get the max values
      MPI_Allreduce(ts.data(), max.data(), ts.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      // Get the min values
      MPI_Allreduce(ts.data(), min.data(), ts.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

      // Get the mean values
      MPI_Allreduce(MPI_IN_PLACE, ts.data(), ts.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      int size;
      // Get the number of CPUs involved
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      // Compute mean times per CPU
      ts /= static_cast<MHDFloat>(size);
   }

}
}
}
