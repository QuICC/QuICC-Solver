/**
 * @file MeanVelocityX.cpp
 * @brief Source of the Mean velocity X physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanVelocityX.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanVelocityX::sTag()
   {
      return "mean_velocityx";
   }

   std::string MeanVelocityX::sFormatted()
   {
      return "Mean velocity X";
   }

   MeanVelocityX::MeanVelocityX()
      : IRegisterId<MeanVelocityX>(MeanVelocityX::sTag(), MeanVelocityX::sFormatted())
   {
   }

}
}
