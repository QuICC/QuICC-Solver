/**
 * @file MeanMagneticZ.cpp
 * @brief Source of the Mean magnetic Z physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanMagneticZ.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanMagneticZ::sTag()
   {
      return "mean_magneticz";
   }

   std::string MeanMagneticZ::sFormatted()
   {
      return "Mean magnetic Z";
   }

   MeanMagneticZ::MeanMagneticZ()
      : IRegisterId<MeanMagneticZ>(MeanMagneticZ::sTag(), MeanMagneticZ::sFormatted())
   {
   }

}
}
