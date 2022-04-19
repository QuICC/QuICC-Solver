/**
 * @file FluctVelocity.cpp
 * @brief Source of the Fluctuating velocity physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctVelocity.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctVelocity::sTag()
   {
      return "fluct_velocity";
   }

   std::string FluctVelocity::sFormatted()
   {
      return "Fluctuating velocity";
   }

   FluctVelocity::FluctVelocity()
      : IRegisterId<FluctVelocity>(FluctVelocity::sTag(), FluctVelocity::sFormatted())
   {
   }

}
}
