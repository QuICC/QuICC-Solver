/**
 * @file MeanMagnetic.cpp
 * @brief Source of the Mean magnetic physical name
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/PhysicalNames/MeanMagnetic.hpp"

// Project includes
//

namespace QuICC {

namespace PhysicalNames {

   std::string MeanMagnetic::sTag()
   {
      return "mean_magnetic";
   }

   std::string MeanMagnetic::sFormatted()
   {
      return "Mean magnetic";
   }

   MeanMagnetic::MeanMagnetic()
      : IRegisterId<MeanMagnetic>(MeanMagnetic::sTag(), MeanMagnetic::sFormatted())
   {
   }

}
}
