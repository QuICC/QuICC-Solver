/**
 * @file DzMeanTemperature.cpp
 * @brief Source of the D_z mean temperature physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/DzMeanTemperature.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string DzMeanTemperature::sTag()
   {
      return "dz_meantemperature";
   }

   std::string DzMeanTemperature::sFormatted()
   {
      return "D_z mean temperature";
   }

   DzMeanTemperature::DzMeanTemperature()
      : IRegisterId<DzMeanTemperature>(DzMeanTemperature::sTag(), DzMeanTemperature::sFormatted())
   {
   }

}
}
