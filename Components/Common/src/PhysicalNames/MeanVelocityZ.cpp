/**
 * @file MeanVelocityZ.cpp
 * @brief Source of the Mean velocity Z physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanVelocityZ.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanVelocityZ::sTag()
   {
      return "mean_velocityz";
   }

   std::string MeanVelocityZ::sFormatted()
   {
      return "Mean velocity Z";
   }

   MeanVelocityZ::MeanVelocityZ()
      : IRegisterId<MeanVelocityZ>(MeanVelocityZ::sTag(), MeanVelocityZ::sFormatted())
   {
   }

}
}
