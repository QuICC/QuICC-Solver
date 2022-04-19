/**
 * @file DxMeanTemperature.cpp
 * @brief Source of the D_x mean temperature physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/DxMeanTemperature.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string DxMeanTemperature::sTag()
   {
      return "dx_meantemperature";
   }

   std::string DxMeanTemperature::sFormatted()
   {
      return "D_x mean temperature";
   }

   DxMeanTemperature::DxMeanTemperature()
      : IRegisterId<DxMeanTemperature>(DxMeanTemperature::sTag(), DxMeanTemperature::sFormatted())
   {
   }

}
}
