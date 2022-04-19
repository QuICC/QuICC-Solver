/**
 * @file Eady.cpp
 * @brief Source of the Eady nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Eady.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Eady::sTag()
   {
      return "eady";
   }

   std::string Eady::sFormatted()
   {
      return "Eady";
   }

   Eady::Eady(const MHDFloat value)
      : IRegisterId<Eady>(value, Eady::sTag(), Eady::sFormatted())
   {
   }

}
}
