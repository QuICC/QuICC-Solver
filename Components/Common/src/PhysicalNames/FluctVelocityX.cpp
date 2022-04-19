/**
 * @file FluctVelocityX.cpp
 * @brief Source of the Fluctuating velocity X physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctVelocityX.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctVelocityX::sTag()
   {
      return "fluct_velocityx";
   }

   std::string FluctVelocityX::sFormatted()
   {
      return "Fluctuating velocity X";
   }

   FluctVelocityX::FluctVelocityX()
      : IRegisterId<FluctVelocityX>(FluctVelocityX::sTag(), FluctVelocityX::sFormatted())
   {
   }

}
}
