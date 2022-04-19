/**
 * @file Sigma.cpp
 * @brief Source of the Sigma nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Sigma.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Sigma::sTag()
   {
      return "sigma";
   }

   std::string Sigma::sFormatted()
   {
      return "Sigma";
   }

   Sigma::Sigma(const MHDFloat value)
      : IRegisterId<Sigma>(value, Sigma::sTag(), Sigma::sFormatted())
   {
   }

}
}
