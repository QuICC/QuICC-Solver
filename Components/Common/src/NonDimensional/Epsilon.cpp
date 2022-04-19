/**
 * @file Epsilon.cpp
 * @brief Source of the Epsilon nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Epsilon.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Epsilon::sTag()
   {
      return "epsilon";
   }

   std::string Epsilon::sFormatted()
   {
      return "Epsilon";
   }

   Epsilon::Epsilon(const MHDFloat value)
      : IRegisterId<Epsilon>(value, Epsilon::sTag(), Epsilon::sFormatted())
   {
   }

}
}
