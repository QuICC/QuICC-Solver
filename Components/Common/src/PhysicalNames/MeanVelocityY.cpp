/**
 * @file MeanVelocityY.cpp
 * @brief Source of the Mean velocity Y physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanVelocityY.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanVelocityY::sTag()
   {
      return "mean_velocityy";
   }

   std::string MeanVelocityY::sFormatted()
   {
      return "Mean velocity Y";
   }

   MeanVelocityY::MeanVelocityY()
      : IRegisterId<MeanVelocityY>(MeanVelocityY::sTag(), MeanVelocityY::sFormatted())
   {
   }

}
}
