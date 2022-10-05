/**
 * @file RadialPower.cpp
 * @brief Source of the reductor transform operator Reductor::RadialPower
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/RadialPower.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string RadialPower::sTag()
   {
      return "Red::RadialPower";
   }

   std::string RadialPower::sFormatted()
   {
      return "Reductor::RadialPower";
   }

   RadialPower::RadialPower()
      : IRegisterId<RadialPower>(RadialPower::sTag(), RadialPower::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
