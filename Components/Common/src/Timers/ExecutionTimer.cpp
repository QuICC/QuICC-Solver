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
      : TimerMacro(autostart), mTimes(Array::Zero(TOTAL + 1))
   {
   }

   ExecutionTimer::~ExecutionTimer()
   {
   }

   void ExecutionTimer::update(BreakPoint point)
   {
      // Increment measured time
      this->mTimes(point) += this->time();

      // Unless TOTAL is request, add time to total execution time
      if(point != TOTAL)
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
      oss << "Initialisation [s]: " << std::fixed << std::setprecision(digits) << this->mTimes(INIT);
      if(max.size() != 0)
      {
         oss << " / " << max(INIT) << " / " << min(INIT);
      }
      Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      // Output prerun time
      oss << "PreRun [s]: " << std::fixed << std::setprecision(digits) << this->mTimes(PRERUN);
      if(max.size() != 0)
      {
         oss << " / " << max(PRERUN) << " / " << min(PRERUN);
      }
      Tools::Formatter::printCentered(stream, oss.str(),' ', base);
      oss.str("");

      // Output computation time
      oss << "Computation [s]: " << std::fixed << std::setprecision(digits) << this->mTimes(RUN);
      if(max.size() != 0)
      {
         oss << " / " << max(RUN) << " / " << min(RUN);
      }
      Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      // Output postrun time
      oss << "PostRun [s]: " << std::fixed << std::setprecision(digits) << this->mTimes(POSTRUN);
      if(max.size() != 0)
      {
         oss << " / " << max(POSTRUN) << " / " << min(POSTRUN);
      }
      Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      Tools::Formatter::printLine(stream, '-');

      // Output total execution time
      oss << "Total execution time: " << std::fixed << std::setprecision(digits) << this->mTimes(TOTAL);
      if(max.size() != 0)
      {
         oss << " / " << max(TOTAL) << " / " << min(TOTAL);
      }
      oss << " seconds";
      Tools::Formatter::printCentered(stream, oss.str(), ' ', base);
      oss.str("");

      Tools::Formatter::printLine(stream, '*');
   }

}
