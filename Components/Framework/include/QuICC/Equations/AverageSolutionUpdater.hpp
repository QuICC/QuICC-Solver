/**
 * @file AverageSolutionUpdater.hpp
 * @brief Solution updater storing averaged value
 */

#pragma once

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "QuICC/Equations/SolutionUpdater.hpp"

namespace QuICC {

namespace Equations {

// THIS IMPLEMENTATION IS NOT COMPLETE. ITS ONLY A PLACEHOLDER

/**
 * @brief Solution updater storing averaged valued
 */
class AverageSolutionUpdater: public SolutionUpdater
{
   public:
      /**
       * @brief ctor
       */
      AverageSolutionUpdater() = default;

      /**
       * @brief dtor
       */
      virtual ~AverageSolutionUpdater() = default;

      /**
       * @brief Simple passthrough
       */
      MHDVariant operator()(const MHDVariant newData, const int, const int, const int) const;
};

MHDVariant AverageSolutionUpdater::operator()(const MHDVariant newData, const int, const int, const int) const
{
   // Only update mean on full timestep
   if(this->mTimeFinished)
   {
      T val = incrementTimeAverage(this->mTimeAvg->point(i, j, k), newData, this->time(), this->mTimestep);
      this->mTimeAvg->setPoint(val, i, j, k);
      return val;
   } else
   {
      return noupdateTimeAverage(this->mTimeAvg->point(i,j,k), newData);
   }
}

} // namespace Equations
} // namespace QuICC
