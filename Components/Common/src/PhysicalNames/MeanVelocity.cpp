/**
 * @file MeanVelocity.cpp
 * @brief Source of the Mean velocity physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanVelocity.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanVelocity::sTag()
   {
      return "mean_velocity";
   }

   std::string MeanVelocity::sFormatted()
   {
      return "Mean velocity";
   }

   MeanVelocity::MeanVelocity()
      : IRegisterId<MeanVelocity>(MeanVelocity::sTag(), MeanVelocity::sFormatted())
   {
   }

}
}
