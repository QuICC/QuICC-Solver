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

   int SpectralTriangularAllL::truncationBwd(const int nN, const int j, const int k)
   {
      return std::max(nN - j/2, MIN_TRUNCATION);
   }

   int SpectralTriangularAllL::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool SpectralTriangularAllL::isOptimal(const int nN, const int maxL)
   {
      return (MIN_TRUNCATION == nN - maxL/2);
   }

} // Tools
} // SpatialScheme
} // QuICC
