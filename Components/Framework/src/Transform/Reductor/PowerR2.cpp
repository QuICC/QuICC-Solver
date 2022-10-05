/**
 * @file PowerR2.cpp
 * @brief Source of the reductor transform operator Reductor::PowerR2
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/PowerR2.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string PowerR2::sTag()
   {
      return "Red::PowerR2";
   }

   std::string PowerR2::sFormatted()
   {
      return "Reductor::PowerR2";
   }

   PowerR2::PowerR2()
      : IRegisterId<PowerR2>(PowerR2::sTag(), PowerR2::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
