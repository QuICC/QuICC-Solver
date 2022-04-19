/**
 * @file VelocityS.cpp
 * @brief Source of the Velocity S physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/VelocityS.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string VelocityS::sTag()
   {
      return "velocitys";
   }

   std::string VelocityS::sFormatted()
   {
      return "Velocity S";
   }

   VelocityS::VelocityS()
      : IRegisterId<VelocityS>(VelocityS::sTag(), VelocityS::sFormatted())
   {
   }

}
}
