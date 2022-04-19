/**
 * @file VelocityZ.cpp
 * @brief Source of the Velocity Z physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/VelocityZ.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string VelocityZ::sTag()
   {
      return "velocityz";
   }

   std::string VelocityZ::sFormatted()
   {
      return "Velocity Z";
   }

   VelocityZ::VelocityZ()
      : IRegisterId<VelocityZ>(VelocityZ::sTag(), VelocityZ::sFormatted())
   {
   }

}
}
