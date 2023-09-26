/** 
 * @file SpectralTrapezoidalSH.cpp
 * @brief Source of the tools for trapezoidal + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/SpectralTrapezoidalSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   SpectralTrapezoidalSH::SpectralTrapezoidalSH()
      : IBaseSH(static_cast<int>(MinimalTruncation::Triangular))
   {
   }

   SpectralTrapezoidalSH::SpectralTrapezoidalSH(const int min)
      : IBaseSH(min)
   {
   }

   int SpectralTrapezoidalSH::truncationFwd(const int nN, const int j, const int k)
   {
      return this->truncationBwd(nN, j, k);
   }

   int SpectralTrapezoidalSH::truncationBwd(const int nN, const int j, const int k)
   {
      return std::max(nN - k/2, this->min());
   }

   int SpectralTrapezoidalSH::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool SpectralTrapezoidalSH::isOptimal(const int nN, const int maxL)
   {
      return (nN - maxL/2 >= this->min());
   }

} // Tools
} // SpatialScheme
} // QuICC
