/**
 * @file ModElsasser.cpp
 * @brief Source of the modified Elsasser nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/ModElsasser.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string ModElsasser::sTag()
   {
      return "modified_elsasser";
   }

   std::string ModElsasser::sFormatted()
   {
      return "modified Elsasser";
   }

   ModElsasser::ModElsasser(const MHDFloat value)
      : IRegisterId<ModElsasser>(value, ModElsasser::sTag(), ModElsasser::sFormatted())
   {
   }

}
}
