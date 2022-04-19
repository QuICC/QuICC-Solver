/**
 * @file Vorticity.cpp
 * @brief Source of the Vorticity physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Vorticity.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Vorticity::sTag()
   {
      return "vorticity";
   }

   std::string Vorticity::sFormatted()
   {
      return "Vorticity";
   }

   Vorticity::Vorticity()
      : IRegisterId<Vorticity>(Vorticity::sTag(), Vorticity::sFormatted())
   {
   }

}
}
