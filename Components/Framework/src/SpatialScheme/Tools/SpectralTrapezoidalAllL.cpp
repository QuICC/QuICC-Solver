/** 
 * @file SpectralTrapezoidalAllL.cpp
 * @brief Source of the tools for trapezoidal + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/SpectralTrapezoidalAllL.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int SpectralTrapezoidalAllL::truncationBwd(const int nN, const int j, const int k)
   {
      return std::max(nN - j/2, MIN_TRUNCATION);
   }

   int SpectralTrapezoidalAllL::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool SpectralTrapezoidalAllL::isOptimal(const int nN, const int maxL)
   {
      return (this->truncationBwd(nN, maxL, 0) > MIN_TRUNCATION);
   }

} // Tools
} // SpatialScheme
} // QuICC
