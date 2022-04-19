/**
 * @file Upsilon.cpp
 * @brief Source of the Upsilon nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Upsilon.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Upsilon::sTag()
   {
      return "upsilon";
   }

   std::string Upsilon::sFormatted()
   {
      return "Upsilon";
   }

   Upsilon::Upsilon(const MHDFloat value)
      : IRegisterId<Upsilon>(value, Upsilon::sTag(), Upsilon::sFormatted())
   {
   }

}
}
