/**
 * @file Eta.cpp
 * @brief Source of the Eta nondimensional number
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/NonDimensional/Eta.hpp"

// Project includes
//

namespace QuICC {

namespace NonDimensional {

   std::string Eta::sTag()
   {
      return "eta";
   }

   std::string Eta::sFormatted()
   {
      return "Eta";
   }

   Eta::Eta(const MHDFloat value)
      : IRegisterId<Eta>(value, Eta::sTag(), Eta::sFormatted())
   {
   }

}
}
