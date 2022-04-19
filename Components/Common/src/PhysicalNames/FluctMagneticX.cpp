/**
 * @file FluctMagneticX.cpp
 * @brief Source of the Fluctuating magnetic X physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctMagneticX.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctMagneticX::sTag()
   {
      return "fluct_magneticx";
   }

   std::string FluctMagneticX::sFormatted()
   {
      return "Fluctuating magnetic X";
   }

   FluctMagneticX::FluctMagneticX()
      : IRegisterId<FluctMagneticX>(FluctMagneticX::sTag(), FluctMagneticX::sFormatted())
   {
   }

}
}
