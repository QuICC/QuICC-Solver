/**
 * @file Zeta.cpp
 * @brief Source of the Zeta nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Zeta.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Zeta::sTag()
   {
      return "zeta";
   }

   std::string Zeta::sFormatted()
   {
      return "Zeta";
   }

   Zeta::Zeta(const MHDFloat value)
      : IRegisterId<Zeta>(value, Zeta::sTag(), Zeta::sFormatted())
   {
   }

}
}
