/**
 * @file FluctVelocityY.cpp
 * @brief Source of the Fluctuating velocity Y physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctVelocityY.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctVelocityY::sTag()
   {
      return "fluct_velocityy";
   }

   std::string FluctVelocityY::sFormatted()
   {
      return "Fluctuating velocity Y";
   }

   FluctVelocityY::FluctVelocityY()
      : IRegisterId<FluctVelocityY>(FluctVelocityY::sTag(), FluctVelocityY::sFormatted())
   {
   }

}
}
