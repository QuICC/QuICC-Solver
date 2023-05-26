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

   const int SpectralTriangularSH::MIN_TRUNCATION = 3;

   int SpectralTriangularSH::truncationFwd(const int nN, const int l)
   {
      return this->truncationBwd(nN, l);
   }

   int SpectralTriangularSH::truncationBwd(const int nN, const int l)
   {
      return std::max(nN - l/2, MIN_TRUNCATION);
   }

   int SpectralTriangularSH::index(const int i, const int l)
   {
      return i;
   }

   bool SpectralTriangularSH::isOptimal(const int nN, const int maxL)
   {
      return ((nN - maxL/2) == this->truncationBwd(nN, maxL));
   }

} // Tools
} // SpatialScheme
} // QuICC
