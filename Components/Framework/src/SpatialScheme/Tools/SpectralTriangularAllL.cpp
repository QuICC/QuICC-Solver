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

   const int SpectralTriangularAllL::MIN_TRUNCATION = 3;

   int SpectralTriangularAllL::truncationBwd(const int nN, const int l)
   {
      return std::max(nN - l/2, MIN_TRUNCATION);
   }

   int SpectralTriangularAllL::index(const int i, const int k)
   {
      return i;
   }

   bool SpectralTriangularAllL::isOptimal(const int nN, const int maxL)
   {
      return ((nN - maxL/2) == this->truncationBwd(nN, maxL));
   }

} // Tools
} // SpatialScheme
} // QuICC
