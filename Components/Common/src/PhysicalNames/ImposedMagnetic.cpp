/**
 * @file ImposedMagnetic.cpp
 * @brief Source of the ImposedMagnetic physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/ImposedMagnetic.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string ImposedMagnetic::sTag()
   {
      return "imposed_magnetic";
   }

   std::string ImposedMagnetic::sFormatted()
   {
      return "Imposed Magnetic";
   }

   ImposedMagnetic::ImposedMagnetic()
      : IRegisterId<ImposedMagnetic>(ImposedMagnetic::sTag(), ImposedMagnetic::sFormatted())
   {
   }

}
}
