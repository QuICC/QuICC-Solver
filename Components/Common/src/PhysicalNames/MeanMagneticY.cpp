/**
 * @file MeanMagneticY.cpp
 * @brief Source of the Mean magnetic Y physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanMagneticY.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanMagneticY::sTag()
   {
      return "mean_magneticy";
   }

   std::string MeanMagneticY::sFormatted()
   {
      return "Mean magnetic Y";
   }

   MeanMagneticY::MeanMagneticY()
      : IRegisterId<MeanMagneticY>(MeanMagneticY::sTag(), MeanMagneticY::sFormatted())
   {
   }

}
}
