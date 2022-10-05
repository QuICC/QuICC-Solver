/**
 * @file Power.cpp
 * @brief Source of the reductor transform operator Reductor::Power
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/Power.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string Power::sTag()
   {
      return "Red::Power";
   }

   std::string Power::sFormatted()
   {
      return "Reductor::Power";
   }

   Power::Power()
      : IRegisterId<Power>(Power::sTag(), Power::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
