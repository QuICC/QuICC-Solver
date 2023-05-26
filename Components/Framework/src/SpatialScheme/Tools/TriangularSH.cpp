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

   const int TriangularSH::MIN_TRUNCATION = 3;

   int TriangularSH::truncationFwd(const int nN, const int l)
   {
      return nN;
   }

   int TriangularSH::truncationBwd(const int nN, const int l)
   {
      return std::max(nN - l/2, MIN_TRUNCATION);
   }

   int TriangularSH::index(const int i, const int l)
   {
      return i;
   }

   bool TriangularSH::isOptimal(const int nN, const int maxL)
   {
      return ((nN - maxL/2) == this->truncationBwd(nN, maxL));
   }

} // Tools
} // SpatialScheme
} // QuICC
