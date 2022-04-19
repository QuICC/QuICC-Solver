/**
 * @file VelocityX.cpp
 * @brief Source of the Velocity X physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/VelocityX.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string VelocityX::sTag()
   {
      return "velocityx";
   }

   std::string VelocityX::sFormatted()
   {
      return "Velocity X";
   }

   VelocityX::VelocityX()
      : IRegisterId<VelocityX>(VelocityX::sTag(), VelocityX::sFormatted())
   {
   }

}
}
