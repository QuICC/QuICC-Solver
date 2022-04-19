/**
 * @file Iota.cpp
 * @brief Source of the Iota nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Iota.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Iota::sTag()
   {
      return "iota";
   }

   std::string Iota::sFormatted()
   {
      return "Iota";
   }

   Iota::Iota(const MHDFloat value)
      : IRegisterId<Iota>(value, Iota::sTag(), Iota::sFormatted())
   {
   }

}
}
