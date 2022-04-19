/**
 * @file Psi.cpp
 * @brief Source of the Psi nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Psi.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Psi::sTag()
   {
      return "psi";
   }

   std::string Psi::sFormatted()
   {
      return "Psi";
   }

   Psi::Psi(const MHDFloat value)
      : IRegisterId<Psi>(value, Psi::sTag(), Psi::sFormatted())
   {
   }

}
}
