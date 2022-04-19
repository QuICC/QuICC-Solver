/**
 * @file Alpha.cpp
 * @brief Source of the Alpha nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Alpha.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Alpha::sTag()
   {
      return "alpha";
   }

   std::string Alpha::sFormatted()
   {
      return "Alpha";
   }

   Alpha::Alpha(const MHDFloat value)
      : IRegisterId<Alpha>(value, Alpha::sTag(), Alpha::sFormatted())
   {
   }

}
}
