/**
 * @file VorticityY.cpp
 * @brief Source of the Vorticity Y physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/VorticityY.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string VorticityY::sTag()
   {
      return "vorticityy";
   }

   std::string VorticityY::sFormatted()
   {
      return "Vorticity Y";
   }

   VorticityY::VorticityY()
      : IRegisterId<VorticityY>(VorticityY::sTag(), VorticityY::sFormatted())
   {
   }

}
}
