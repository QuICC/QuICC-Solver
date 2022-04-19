/**
 * @file Tau.cpp
 * @brief Source of the Tau nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Tau.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Tau::sTag()
   {
      return "tau";
   }

   std::string Tau::sFormatted()
   {
      return "Tau";
   }

   Tau::Tau(const MHDFloat value)
      : IRegisterId<Tau>(value, Tau::sTag(), Tau::sFormatted())
   {
   }

}
}
