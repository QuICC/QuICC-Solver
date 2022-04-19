/**
 * @file Heating.cpp
 * @brief Source of the Heating nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Heating.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Heating::sTag()
   {
      return "heating";
   }

   std::string Heating::sFormatted()
   {
      return "Heating";
   }

   Heating::Heating(const MHDFloat value)
      : IRegisterId<Heating>(value, Heating::sTag(), Heating::sFormatted())
   {
   }

}
}
