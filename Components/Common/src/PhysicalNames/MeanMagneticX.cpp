/**
 * @file MeanMagneticX.cpp
 * @brief Source of the Mean magnetic X physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanMagneticX.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanMagneticX::sTag()
   {
      return "mean_magneticx";
   }

   std::string MeanMagneticX::sFormatted()
   {
      return "Mean magnetic X";
   }

   MeanMagneticX::MeanMagneticX()
      : IRegisterId<MeanMagneticX>(MeanMagneticX::sTag(), MeanMagneticX::sFormatted())
   {
   }

}
}
