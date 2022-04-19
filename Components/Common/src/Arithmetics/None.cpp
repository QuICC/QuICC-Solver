/**
 * @file None.cpp
 * @brief Source of the None Arithmetics
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Arithmetics/None.hpp"

// Project includes
//

namespace QuICC {

namespace Arithmetics {

   std::string None::sTag()
   {
      return "none";
   }

   std::string None::sFormatted()
   {
      return "None";
   }

   None::None()
      : IRegisterId<None>(None::sTag(), None::sFormatted())
   {
   }

}
}
