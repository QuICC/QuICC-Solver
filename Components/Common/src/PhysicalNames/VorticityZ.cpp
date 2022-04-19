/**
 * @file VorticityZ.cpp
 * @brief Source of the Vorticity Z physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/VorticityZ.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string VorticityZ::sTag()
   {
      return "vorticityz";
   }

   std::string VorticityZ::sFormatted()
   {
      return "Vorticity Z";
   }

   VorticityZ::VorticityZ()
      : IRegisterId<VorticityZ>(VorticityZ::sTag(), VorticityZ::sFormatted())
   {
   }

}
}
