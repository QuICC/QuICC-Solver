/**
 * @file Kappa.cpp
 * @brief Source of the Kappa nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Kappa.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Kappa::sTag()
   {
      return "kappa";
   }

   std::string Kappa::sFormatted()
   {
      return "Kappa";
   }

   Kappa::Kappa(const MHDFloat value)
      : IRegisterId<Kappa>(value, Kappa::sTag(), Kappa::sFormatted())
   {
   }

}
}
