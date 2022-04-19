/**
 * @file VelocityY.cpp
 * @brief Source of the Velocity Y physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/VelocityY.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string VelocityY::sTag()
   {
      return "velocityy";
   }

   std::string VelocityY::sFormatted()
   {
      return "Velocity Y";
   }

   VelocityY::VelocityY()
      : IRegisterId<VelocityY>(VelocityY::sTag(), VelocityY::sFormatted())
   {
   }

}
}
