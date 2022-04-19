/**
 * @file Rescaled.cpp
 * @brief Source of the Rescaled nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Rescaled.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Rescaled::sTag()
   {
      return "rescaled";
   }

   std::string Rescaled::sFormatted()
   {
      return "Rescaled";
   }

   Rescaled::Rescaled(const MHDFloat value)
      : IRegisterId<Rescaled>(value, Rescaled::sTag(), Rescaled::sFormatted())
   {
   }

}
}
