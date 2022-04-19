/**
 * @file Theta.cpp
 * @brief Source of the Theta nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Theta.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Theta::sTag()
   {
      return "theta";
   }

   std::string Theta::sFormatted()
   {
      return "Theta";
   }

   Theta::Theta(const MHDFloat value)
      : IRegisterId<Theta>(value, Theta::sTag(), Theta::sFormatted())
   {
   }

}
}
