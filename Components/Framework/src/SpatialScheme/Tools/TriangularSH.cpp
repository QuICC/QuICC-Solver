/** 
 * @file TriangularSH.cpp
 * @brief Source of the tools for triangular + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/TriangularSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int TriangularSH::truncationFwd(const int nN, const int j, const int k)
   {
      return nN;
   }

   int TriangularSH::truncationBwd(const int nN, const int j, const int k)
   {
      return std::max(nN - k/2, MIN_TRUNCATION);
   }

   int TriangularSH::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool TriangularSH::isOptimal(const int nN, const int maxL)
   {
      return (MIN_TRUNCATION == nN - maxL/2);
   }

} // Tools
} // SpatialScheme
} // QuICC
