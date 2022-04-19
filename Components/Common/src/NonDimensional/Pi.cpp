/**
 * @file Pi.cpp
 * @brief Source of the Pi nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Pi.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Pi::sTag()
   {
      return "pi";
   }

   std::string Pi::sFormatted()
   {
      return "Pi";
   }

   Pi::Pi(const MHDFloat value)
      : IRegisterId<Pi>(value, Pi::sTag(), Pi::sFormatted())
   {
   }

}
}
