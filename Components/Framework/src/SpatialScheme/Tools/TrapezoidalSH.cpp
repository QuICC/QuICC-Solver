/** 
 * @file TrapezoidalSH.cpp
 * @brief Source of the tools for trapezoidal + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/TrapezoidalSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int TrapezoidalSH::truncationFwd(const int nN, const int j, const int k)
   {
      return nN;
   }

   int TrapezoidalSH::truncationBwd(const int nN, const int j, const int k)
   {
      return std::max(nN - k/2, MIN_TRUNCATION);
   }

   int TrapezoidalSH::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool TrapezoidalSH::isOptimal(const int nN, const int maxL)
   {
      return (this->truncationBwd(nN, 0, maxL) > MIN_TRUNCATION);
   }

} // Tools
} // SpatialScheme
} // QuICC
