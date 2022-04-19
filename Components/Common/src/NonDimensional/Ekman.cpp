/**
 * @file Ekman.cpp
 * @brief Source of the Ekman nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Ekman.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Ekman::sTag()
   {
      return "ekman";
   }

   std::string Ekman::sFormatted()
   {
      return "Ekman";
   }

   Ekman::Ekman(const MHDFloat value)
      : IRegisterId<Ekman>(value, Ekman::sTag(), Ekman::sFormatted())
   {
   }

}
}
