/**
 * @file Velocity.cpp
 * @brief Source of the Velocity physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/Velocity.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string Velocity::sTag()
   {
      return "velocity";
   }

   std::string Velocity::sFormatted()
   {
      return "Velocity";
   }

   Velocity::Velocity()
      : IRegisterId<Velocity>(Velocity::sTag(), Velocity::sFormatted())
   {
   }

}
}
