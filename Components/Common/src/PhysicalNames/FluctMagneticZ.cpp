/**
 * @file FluctMagneticZ.cpp
 * @brief Source of the Fluctuating magnetic Z physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctMagneticZ.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctMagneticZ::sTag()
   {
      return "fluct_magneticz";
   }

   std::string FluctMagneticZ::sFormatted()
   {
      return "Fluctuating magnetic Z";
   }

   FluctMagneticZ::FluctMagneticZ()
      : IRegisterId<FluctMagneticZ>(FluctMagneticZ::sTag(), FluctMagneticZ::sFormatted())
   {
   }

}
}
