/**
 * @file PowerD1.cpp
 * @brief Source of the reductor transform operator Reductor::PowerD1
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/PowerD1.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string PowerD1::sTag()
   {
      return "Red::PowerD1";
   }

   std::string PowerD1::sFormatted()
   {
      return "Reductor::PowerD1";
   }

   PowerD1::PowerD1()
      : IRegisterId<PowerD1>(PowerD1::sTag(), PowerD1::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
