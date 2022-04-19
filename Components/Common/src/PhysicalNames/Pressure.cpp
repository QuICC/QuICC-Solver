/**
 * @file Pressure.cpp
 * @brief Source of the Pressure physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Pressure.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Pressure::sTag()
   {
      return "pressure";
   }

   std::string Pressure::sFormatted()
   {
      return "Pressure";
   }

   Pressure::Pressure()
      : IRegisterId<Pressure>(Pressure::sTag(), Pressure::sFormatted())
   {
   }

}
}
