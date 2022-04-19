/**
 * @file VorticityX.cpp
 * @brief Source of the Vorticity X physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/VorticityX.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string VorticityX::sTag()
   {
      return "vorticityx";
   }

   std::string VorticityX::sFormatted()
   {
      return "Vorticity X";
   }

   VorticityX::VorticityX()
      : IRegisterId<VorticityX>(VorticityX::sTag(), VorticityX::sFormatted())
   {
   }

}
}
