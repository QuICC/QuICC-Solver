/**
 * @file Phi.cpp
 * @brief Source of the Phi nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Phi.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Phi::sTag()
   {
      return "phi";
   }

   std::string Phi::sFormatted()
   {
      return "Phi";
   }

   Phi::Phi(const MHDFloat value)
      : IRegisterId<Phi>(value, Phi::sTag(), Phi::sFormatted())
   {
   }

}
}
