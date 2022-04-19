/**
 * @file Before.cpp
 * @brief Source of the Before SolveTiming
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SolveTiming/Before.hpp"

// Project includes
//

namespace QuICC {

namespace SolveTiming {

   std::string Before::sTag()
   {
      return "before";
   }

   std::string Before::sFormatted()
   {
      return "Before";
   }

   Before::Before()
      : IRegisterId<Before>(Before::sTag(), Before::sFormatted())
   {
   }

}
}
