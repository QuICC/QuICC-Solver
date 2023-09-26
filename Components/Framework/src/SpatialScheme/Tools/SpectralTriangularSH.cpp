/** 
 * @file SpectralTriangularSH.cpp
 * @brief Source of the tools for triangular + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/SpectralTriangularSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   SpectralTriangularSH::SpectralTriangularSH(const int min)
      : SpectralTrapezoidalSH(min)
   {
   }

   bool SpectralTriangularSH::isOptimal(const int nN, const int maxL)
   {
      return (this->min() == nN - maxL/2);
   }

} // Tools
} // SpatialScheme
} // QuICC
