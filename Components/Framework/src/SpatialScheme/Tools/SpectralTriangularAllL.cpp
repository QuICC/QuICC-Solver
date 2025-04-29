/** 
 * @file SpectralTriangularAllL.cpp
 * @brief Source of the tools for triangular + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/SpectralTriangularAllL.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   SpectralTriangularAllL::SpectralTriangularAllL(const int min)
      : SpectralTrapezoidalAllL(min)
   {
   }

   bool SpectralTriangularAllL::isOptimal(const int nN, const int maxL)
   {
      return (this->min() == nN - maxL/2);
   }

} // Tools
} // SpatialScheme
} // QuICC
