/**
 * @file VorticityS.cpp
 * @brief Source of the Vorticity S physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/VorticityS.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string VorticityS::sTag()
   {
      return "vorticitys";
   }

   std::string VorticityS::sFormatted()
   {
      return "Vorticity S";
   }

   VorticityS::VorticityS()
      : IRegisterId<VorticityS>(VorticityS::sTag(), VorticityS::sFormatted())
   {
   }

}
}
