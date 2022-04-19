/**
 * @file MagneticX.cpp
 * @brief Source of the Magnetic X physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MagneticX.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MagneticX::sTag()
   {
      return "magneticx";
   }

   std::string MagneticX::sFormatted()
   {
      return "Magnetic X";
   }

   MagneticX::MagneticX()
      : IRegisterId<MagneticX>(MagneticX::sTag(), MagneticX::sFormatted())
   {
   }

}
}
