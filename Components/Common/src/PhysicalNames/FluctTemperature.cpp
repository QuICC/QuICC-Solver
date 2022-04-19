/**
 * @file FluctTemperature.cpp
 * @brief Source of the Fluctuating temperature physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/FluctTemperature.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string FluctTemperature::sTag()
   {
      return "fluct_temperature";
   }

   std::string FluctTemperature::sFormatted()
   {
      return "Fluctuating temperature";
   }

   FluctTemperature::FluctTemperature()
      : IRegisterId<FluctTemperature>(FluctTemperature::sTag(), FluctTemperature::sFormatted())
   {
   }

}
}
