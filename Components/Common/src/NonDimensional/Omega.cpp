/**
 * @file Omega.cpp
 * @brief Source of the Omega nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Omega.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Omega::sTag()
   {
      return "omega";
   }

   std::string Omega::sFormatted()
   {
      return "Omega";
   }

   Omega::Omega(const MHDFloat value)
      : IRegisterId<Omega>(value, Omega::sTag(), Omega::sFormatted())
   {
   }

}
}
