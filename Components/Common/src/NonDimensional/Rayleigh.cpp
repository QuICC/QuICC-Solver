/**
 * @file Rayleigh.cpp
 * @brief Source of the Rayleigh nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Rayleigh.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Rayleigh::sTag()
   {
      return "rayleigh";
   }

   std::string Rayleigh::sFormatted()
   {
      return "Rayleigh";
   }

   Rayleigh::Rayleigh(const MHDFloat value)
      : IRegisterId<Rayleigh>(value, Rayleigh::sTag(), Rayleigh::sFormatted())
   {
   }

}
}
