/**
 * @file Set.cpp
 * @brief Source of the Set Arithmetics
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Arithmetics/Set.hpp"

// Project includes
//

namespace QuICC {

namespace Arithmetics {

   std::string Set::sTag()
   {
      return "set";
   }

   std::string Set::sFormatted()
   {
      return "Set";
   }

   Set::Set()
      : IRegisterId<Set>(Set::sTag(), Set::sFormatted())
   {
   }

}
}
