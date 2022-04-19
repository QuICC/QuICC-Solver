/**
 * @file Roberts.cpp
 * @brief Source of the Roberts nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Roberts.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Roberts::sTag()
   {
      return "roberts";
   }

   std::string Roberts::sFormatted()
   {
      return "Roberts";
   }

   Roberts::Roberts(const MHDFloat value)
      : IRegisterId<Roberts>(value, Roberts::sTag(), Roberts::sFormatted())
   {
   }

}
}
