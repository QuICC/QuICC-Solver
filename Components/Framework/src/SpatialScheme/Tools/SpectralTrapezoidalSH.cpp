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

   const int SpectralTrapezoidalSH::MIN_TRUNCATION = 3;

   int SpectralTrapezoidalSH::truncationFwd(const int nN, const int l)
   {
      return this->truncationBwd(nN, l);
   }

   int SpectralTrapezoidalSH::truncationBwd(const int nN, const int l)
   {
      return std::max(nN - l/2, MIN_TRUNCATION);
   }

   int SpectralTrapezoidalSH::index(const int i, const int l)
   {
      return i;
   }

   bool SpectralTrapezoidalSH::isOptimal(const int nN, const int maxL)
   {
      return (this->truncationBwd(nN, maxL) > MIN_TRUNCATION);
   }

} // Tools
} // SpatialScheme
} // QuICC
