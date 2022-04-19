/**
 * @file Elsasser.cpp
 * @brief Source of the Elsasser nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Elsasser.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Elsasser::sTag()
   {
      return "elsasser";
   }

   std::string Elsasser::sFormatted()
   {
      return "Elsasser";
   }

   Elsasser::Elsasser(const MHDFloat value)
      : IRegisterId<Elsasser>(value, Elsasser::sTag(), Elsasser::sFormatted())
   {
   }

}
}
