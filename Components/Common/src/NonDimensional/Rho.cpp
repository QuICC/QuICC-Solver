/**
 * @file Rho.cpp
 * @brief Source of the Rho nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Rho.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Rho::sTag()
   {
      return "rho";
   }

   std::string Rho::sFormatted()
   {
      return "Rho";
   }

   Rho::Rho(const MHDFloat value)
      : IRegisterId<Rho>(value, Rho::sTag(), Rho::sFormatted())
   {
   }

}
}
