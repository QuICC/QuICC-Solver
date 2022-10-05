/**
 * @file PowerD1R1.cpp
 * @brief Source of the reductor transform operator Reductor::PowerD1R1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/PowerD1R1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string PowerD1R1::sTag()
   {
      return "Red::PowerD1R1";
   }

   std::string PowerD1R1::sFormatted()
   {
      return "Reductor::PowerD1R1";
   }

   PowerD1R1::PowerD1R1()
      : IRegisterId<PowerD1R1>(PowerD1R1::sTag(), PowerD1R1::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
