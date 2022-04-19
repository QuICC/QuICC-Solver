/**
 * @file Sub.cpp
 * @brief Source of the Sub Arithmetics
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Arithmetics/Sub.hpp"

// Project includes
//

namespace QuICC {

namespace Arithmetics {

   std::string Sub::sTag()
   {
      return "sub";
   }

   std::string Sub::sFormatted()
   {
      return "Sub";
   }

   Sub::Sub()
      : IRegisterId<Sub>(Sub::sTag(), Sub::sFormatted())
   {
   }

}
}
