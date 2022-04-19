/**
 * @file MeanTemperature.cpp
 * @brief Source of the Mean temperature physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanTemperature.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanTemperature::sTag()
   {
      return "mean_temperature";
   }

   std::string MeanTemperature::sFormatted()
   {
      return "Mean temperature";
   }

   MeanTemperature::MeanTemperature()
      : IRegisterId<MeanTemperature>(MeanTemperature::sTag(), MeanTemperature::sFormatted())
   {
   }

}
}
