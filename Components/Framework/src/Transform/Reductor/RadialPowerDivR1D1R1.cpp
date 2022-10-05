/**
 * @file RadialPowerDivR1D1R1.cpp
 * @brief Source of the reductor transform operator Reductor::RadialPowerDivR1D1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/RadialPowerDivR1D1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string RadialPowerDivR1D1R1::sTag()
   {
      return "Red::RadialPowerDivR1D1R1";
   }

   std::string RadialPowerDivR1D1R1::sFormatted()
   {
      return "Reductor::RadialPowerDivR1D1R1";
   }

   RadialPowerDivR1D1R1::RadialPowerDivR1D1R1()
      : IRegisterId<RadialPowerDivR1D1R1>(RadialPowerDivR1D1R1::sTag(), RadialPowerDivR1D1R1::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
