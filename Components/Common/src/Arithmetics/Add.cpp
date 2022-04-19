/**
 * @file Add.cpp
 * @brief Source of the Add Arithmetics
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Arithmetics/Add.hpp"

// Project includes
//

namespace QuICC {

namespace Arithmetics {

   std::string Add::sTag()
   {
      return "add";
   }

   std::string Add::sFormatted()
   {
      return "Add";
   }

   Add::Add()
      : IRegisterId<Add>(Add::sTag(), Add::sFormatted())
   {
   }

}
}
