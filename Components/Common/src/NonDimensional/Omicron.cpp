/**
 * @file Omicron.cpp
 * @brief Source of the Omicron nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Omicron.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Omicron::sTag()
   {
      return "omicron";
   }

   std::string Omicron::sFormatted()
   {
      return "Omicron";
   }

   Omicron::Omicron(const MHDFloat value)
      : IRegisterId<Omicron>(value, Omicron::sTag(), Omicron::sFormatted())
   {
   }

}
}
