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

   TrapezoidalSH::TrapezoidalSH()
      : IBaseSH(static_cast<int>(MinimalTruncation::Triangular))
   {
   }

   TrapezoidalSH::TrapezoidalSH(const int min)
      : IBaseSH(min)
   {
   }

   int TrapezoidalSH::truncationFwd(const int nN, const int j, const int k)
   {
      return nN;
   }

   int TrapezoidalSH::truncationBwd(const int nN, const int j, const int k)
   {
      return std::max(nN - k/2, this->min());
   }

   int TrapezoidalSH::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool TrapezoidalSH::isOptimal(const int nN, const int maxL)
   {
      return (nN - maxL/2 >= this->min());
   }

} // Tools
} // SpatialScheme
} // QuICC
