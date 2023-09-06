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

   const int TrapezoidalSH::MIN_TRUNCATION = 3;

   int TrapezoidalSH::truncationFwd(const int nN, const int l)
   {
      return nN;
   }

   int TrapezoidalSH::truncationBwd(const int nN, const int l)
   {
      return std::max(nN - l/2, MIN_TRUNCATION);
   }

   int TrapezoidalSH::index(const int i, const int l)
   {
      return i;
   }

   bool TrapezoidalSH::isOptimal(const int nN, const int maxL)
   {
      return (this->truncationBwd(nN, maxL) > MIN_TRUNCATION);
   }

} // Tools
} // SpatialScheme
} // QuICC
