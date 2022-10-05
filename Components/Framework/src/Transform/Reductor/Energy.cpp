/**
 * @file Energy.cpp
 * @brief Source of the reductor transform operator Reductor::Energy
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/Energy.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

   std::string Energy::sTag()
   {
      return "Red::Energy";
   }

   std::string Energy::sFormatted()
   {
      return "Reductor::Energy";
   }

   Energy::Energy()
      : IRegisterId<Energy>(Energy::sTag(), Energy::sFormatted())
   {
   }

} // Reductor
} // Transform
} // QuICC
