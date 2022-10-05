/**
 * @file RadialPowerDivR1.cpp
 * @brief Source of the reductor transform operator Reductor::RadialPowerDivR1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/RadialPowerDivR1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string RadialPowerDivR1::sTag()
   {
      return "Red::RadialPowerDivR1";
   }

   std::string RadialPowerDivR1::sFormatted()
   {
      return "Reductor::RadialPowerDivR1";
   }

   RadialPowerDivR1::RadialPowerDivR1()
      : IRegisterId<RadialPowerDivR1>(RadialPowerDivR1::sTag(), RadialPowerDivR1::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
