/**
 * @file Gamma.cpp
 * @brief Source of the Gamma nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Gamma.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Gamma::sTag()
   {
      return "gamma";
   }

   std::string Gamma::sFormatted()
   {
      return "Gamma";
   }

   Gamma::Gamma(const MHDFloat value)
      : IRegisterId<Gamma>(value, Gamma::sTag(), Gamma::sFormatted())
   {
   }

}
}
