/**
 * @file PowerSLAPLR2.cpp
 * @brief Source of the reductor transform operator Reductor::PowerSLAPLR2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/PowerSLAPLR2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string PowerSLAPLR2::sTag()
   {
      return "Red::PowerSLAPLR2";
   }

   std::string PowerSLAPLR2::sFormatted()
   {
      return "Reductor::PowerSLAPLR2";
   }

   PowerSLAPLR2::PowerSLAPLR2()
      : IRegisterId<PowerSLAPLR2>(PowerSLAPLR2::sTag(), PowerSLAPLR2::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
