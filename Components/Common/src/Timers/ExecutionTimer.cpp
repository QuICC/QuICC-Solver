/** 
 * @file ExecutionTimer.cpp
 * @brief Source of the implementation of an execution timer
 */

// Configuration includes
//

// System includes
//
#ifdef QUICC_MPI
#include <mpi.h>
#endif // QUICC_MPI

// External includes
//

// Class include
//
#include "QuICC/Timers/ExecutionTimer.hpp"

// Project includes
//
#include "QuICC/Tools/Formatter.hpp"

namespace QuICC {

   ExecutionTimer::ExecutionTimer(const bool autostart)
      : TimerMacro(autostart), mTimes(Array::Zero(NBREAKPOINT)), mCount(0)
   {
   }

   ExecutionTimer::~ExecutionTimer()
   {
   }

   void ExecutionTimer::iteration()
   {
      this->mCount++;

      if(this->mCount == 1)
      {
         // Extract timing of timestep at first step
         auto t = this->mTimes(TIMESTEP);
         this->mTimes(FIRST_TIMESTEP) = t;
         this->mTimes(FIRST_RUN) = t;

         // Extract timing of nonlinear at first step
         t = this->mTimes(NONLINEAR);
         this->mTimes(FIRST_NONLINEAR) = t;
         this->mTimes(FIRST_RUN) += t;
      }
   }

   void ExecutionTimer::update(BreakPoint point)
   {
      // Increment measured time
      this->mTimes(point) += this->time();

      // Unless TOTAL is request, add time to total execution time
      if(point < TOTAL)
      {
         this->mTimes(TOTAL) += this->time();
      }
   }

   double ExecutionTimer::queryTime(BreakPoint point) const
   {
      return this->mTimes(point) + this->queryTime();
   }

   void ExecutionTimer::analyze(Array& min, Array& max)
   {
      // Get the "global" times from MPI code
      #ifdef QUICC_MPI
         // Resize the storage
         min.resize(this->mTimes.size());
         max.resize(this->mTimes.size());

         // Get the max values
         MPI_Allreduce(this->mTimes.data(), max.data(), this->mTimes.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

         // Get the min values
         MPI_Allreduce(this->mTimes.data(), min.data(), this->mTimes.size(), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
         // Get the mean values
         MPI_Allreduce(MPI_IN_PLACE, this->mTimes.data(), this->mTimes.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         // Compute mean times per CPU
         int nCpu;
         MPI_Comm_size(MPI_COMM_WORLD, &nCpu);
         this->mTimes /= static_cast<double>(nCpu);
      #endif // QUICC_MPI
   }

   void ExecutionTimer::printInfo(std::ostream& stream)
   {
      // Analyze the data
      Array min;
      Array max;
      this->analyze(min, max);
      int digits = 3;

      // Create nice looking ouput header
      Tools::Formatter::printNewline(stream);
      Tools::Formatter::printLine(stream, '-');
      Tools::Formatter::printCentered(stream, "Execution time information", '*');
      Tools::Formatter::printLine(stream, '-');

      std::stringstream oss;

      // get a nice base for info
      int base = 20;
      oss << std::fixed << std::setprecision(digits) << this->mTimes.maxCoeff();
      base += oss.str().size() + 1;
      oss.str("");

      // Output initialisation time
      this->formatTiming(stream, "Initialisation", INIT, min, max, 0, digits, base);

      // Output prerun time
      this->formatTiming(stream, "PreRun", PRERUN, min, max, 0, digits, base);

      // Output computation time of first step
      this->formatTiming(stream, "First step", FIRST_RUN, min, max, 0, digits, base);

      // Output timestep time of first step
      this->formatTiming(stream, "Timestep", FIRST_TIMESTEP, min, max, 1, digits, base);

      // Output nonlinear time of first step
      this->formatTiming(stream, "Nonlinear", FIRST_NONLINEAR, min, max, 1, digits, base);

      // Output computation time
      this->formatTiming(stream, "Computation", RUN, min, max, 0, digits, base);

      // Output timestep time
      this->formatTiming(stream, "Timestep", TIMESTEP, min, max, 1, digits, base);

      // Output nonlinear time
      this->formatTiming(stream, "Nonlinear", NONLINEAR, min, max, 1, digits, base);

      // Output postrun time
      this->formatTiming(stream, "PostRun", POSTRUN, min, max, 0, digits, base);

      // Print number of iterations
      oss << "Iterations: " << this->mCount;
      Tools::Formatter::printCentered(stream, oss.str(), ' ', base);

      Tools::Formatter::printLine(stream, '-');

      // Output total execution time
      this->formatTiming(stream, "Total", TOTAL, min, max, 0, digits, base);

      Tools::Formatter::printLine(stream, '*');
   }

   void ExecutionTimer::formatTiming(std::ostream& stream, const std::string name, const BreakPoint pt, const Array& min ,const Array& max, const int spaces, const int digits, const int base)
   {
      std::stringstream oss;

      // Output nonlinear time
      oss << std::string(spaces, '-')  << name << " [s]: " << std::fixed << std::setprecision(digits) << this->mTimes(pt);
      if(max.size() != 0)
      {
         oss << " / " << max(pt) << " / " << min(pt);
      }
      Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");
   }

}
