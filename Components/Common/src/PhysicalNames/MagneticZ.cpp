/**
 * @file MagneticZ.cpp
 * @brief Source of the Magnetic Z physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MagneticZ.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MagneticZ::sTag()
   {
      return "magneticz";
   }

   std::string MagneticZ::sFormatted()
   {
      return "Magnetic Z";
   }

   MagneticZ::MagneticZ()
      : IRegisterId<MagneticZ>(MagneticZ::sTag(), MagneticZ::sFormatted())
   {
   }

}
}
