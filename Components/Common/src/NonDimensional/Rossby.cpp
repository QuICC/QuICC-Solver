/**
 * @file Rossby.cpp
 * @brief Source of the Rossby nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Rossby.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Rossby::sTag()
   {
      return "rossby";
   }

   std::string Rossby::sFormatted()
   {
      return "Rossby";
   }

   Rossby::Rossby(const MHDFloat value)
      : IRegisterId<Rossby>(value, Rossby::sTag(), Rossby::sFormatted())
   {
   }

}
}
