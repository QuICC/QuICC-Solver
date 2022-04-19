/**
 * @file Density.cpp
 * @brief Source of the Density physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Density.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Density::sTag()
   {
      return "density";
   }

   std::string Density::sFormatted()
   {
      return "Density";
   }

   Density::Density()
      : IRegisterId<Density>(Density::sTag(), Density::sFormatted())
   {
   }

}
}
