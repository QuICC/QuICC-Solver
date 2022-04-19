/**
 * @file MagneticY.cpp
 * @brief Source of the Magnetic Y physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MagneticY.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MagneticY::sTag()
   {
      return "magneticy";
   }

   std::string MagneticY::sFormatted()
   {
      return "Magnetic Y";
   }

   MagneticY::MagneticY()
      : IRegisterId<MagneticY>(MagneticY::sTag(), MagneticY::sFormatted())
   {
   }

}
}
