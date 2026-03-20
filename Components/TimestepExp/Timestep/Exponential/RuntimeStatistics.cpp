/**
 * @file RuntimeStatistics.cpp
 * @brief Implementation of runtime statistics
 */

// System includes
//
#include <algorithm>
#include <limits>
#include <iostream>

// Project includes
//
#include "Timestep/Exponential/RuntimeStatistics.hpp"

namespace QuICC {

namespace Timestep {

namespace Exponential {

   RuntimeStatistics::RuntimeStatistics()
      : step(0), krystep(0), reject(0), exps(0), m(-1), mMax(std::numeric_limits<int>::min()), mMin(std::numeric_limits<int>::max()), conv(0), convMax(std::numeric_limits<double>::min()), convMin(std::numeric_limits<double>::max())
   {
   }

   void RuntimeStatistics::reset()
   {
      this->step = 0;
      this->krystep = 0;
      this->reject = 0;
      this->exps = 0;
      this->m = -1;
      this->conv = 0;
   }

   void RuntimeStatistics::update()
   {
      this->mMax = std::max(this->mMax, m);
      this->mMin = std::min(this->mMin, m);
      this->convMax = std::max(this->convMax, conv);
      this->convMin = std::min(this->convMin, conv);
   }

   void RuntimeStatistics::printInfo()
   {
      const int n = 3;
      std::cerr
         << std::string(n, ' ') << "steps: " << this->step
         << std::string(n, ' ') << "krylov steps: " << this->krystep
         << std::string(n, ' ') << "reject: " << this->reject
         << std::string(n, ' ') << "exps: " << this->exps
         << std::string(n, ' ') << "m: " << this->m
         << std::string(n, ' ') << "mMax: " << this->mMax
         << std::string(n, ' ') << "mMin: " << this->mMin
         << std::string(n, ' ') << "conv: " << this->conv
         << std::string(n, ' ') << "convMax: " << this->convMax
         << std::string(n, ' ') << "convMin: " << this->convMin
         << std::endl;
   }

} // namespace Exponential
} // namespace Timestep
} // namespace QuICC
