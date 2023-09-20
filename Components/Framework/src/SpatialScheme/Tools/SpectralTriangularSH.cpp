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

   int SpectralTriangularSH::truncationFwd(const int nN, const int j, const int k)
   {
      return this->truncationBwd(nN, j, k);
   }

   int SpectralTriangularSH::truncationBwd(const int nN, const int j, const int k)
   {
      return std::max(nN - k/2, MIN_TRUNCATION);
   }

   int SpectralTriangularSH::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool SpectralTriangularSH::isOptimal(const int nN, const int maxL)
   {
      return (MIN_TRUNCATION == nN - maxL/2);
   }

} // Tools
} // SpatialScheme
} // QuICC
