/**
 * @file Mu.cpp
 * @brief Source of the Mu nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Mu.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Mu::sTag()
   {
      return "mu";
   }

   std::string Mu::sFormatted()
   {
      return "Mu";
   }

   Mu::Mu(const MHDFloat value)
      : IRegisterId<Mu>(value, Mu::sTag(), Mu::sFormatted())
   {
   }

}
}
