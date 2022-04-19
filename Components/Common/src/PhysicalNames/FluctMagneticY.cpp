/**
 * @file FluctMagneticY.cpp
 * @brief Source of the Fluctuating magnetic Y physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctMagneticY.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctMagneticY::sTag()
   {
      return "fluct_magneticy";
   }

   std::string FluctMagneticY::sFormatted()
   {
      return "Fluctuating magnetic Y";
   }

   FluctMagneticY::FluctMagneticY()
      : IRegisterId<FluctMagneticY>(FluctMagneticY::sTag(), FluctMagneticY::sFormatted())
   {
   }

}
}
