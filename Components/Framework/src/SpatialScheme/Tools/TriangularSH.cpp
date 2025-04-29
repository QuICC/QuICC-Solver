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

   TriangularSH::TriangularSH(const int min)
      : TrapezoidalSH(min)
   {
   }

   bool TriangularSH::isOptimal(const int nN, const int maxL)
   {
      return (this->min() == nN - maxL/2);
   }

} // Tools
} // SpatialScheme
} // QuICC
