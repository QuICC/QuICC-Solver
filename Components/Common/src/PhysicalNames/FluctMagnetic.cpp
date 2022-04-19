/**
 * @file FluctMagnetic.cpp
 * @brief Source of the Fluctuating magnetic physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctMagnetic.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctMagnetic::sTag()
   {
      return "fluct_magnetic";
   }

   std::string FluctMagnetic::sFormatted()
   {
      return "Fluctuating magnetic";
   }

   FluctMagnetic::FluctMagnetic()
      : IRegisterId<FluctMagnetic>(FluctMagnetic::sTag(), FluctMagnetic::sFormatted())
   {
   }

}
}
