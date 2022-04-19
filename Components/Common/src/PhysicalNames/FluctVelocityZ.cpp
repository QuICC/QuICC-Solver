/**
 * @file FluctVelocityZ.cpp
 * @brief Source of the Fluctuating velocity Z physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctVelocityZ.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctVelocityZ::sTag()
   {
      return "fluct_velocityz";
   }

   std::string FluctVelocityZ::sFormatted()
   {
      return "Fluctuating velocity Z";
   }

   FluctVelocityZ::FluctVelocityZ()
      : IRegisterId<FluctVelocityZ>(FluctVelocityZ::sTag(), FluctVelocityZ::sFormatted())
   {
   }

}
}
