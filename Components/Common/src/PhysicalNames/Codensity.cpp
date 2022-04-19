/**
 * @file Codensity.cpp
 * @brief Source of the Codensity physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Codensity.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Codensity::sTag()
   {
      return "codensity";
   }

   std::string Codensity::sFormatted()
   {
      return "Codensity";
   }

   Codensity::Codensity()
      : IRegisterId<Codensity>(Codensity::sTag(), Codensity::sFormatted())
   {
   }

}
}
