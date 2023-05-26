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

   const int SpectralTrapezoidalAllL::MIN_TRUNCATION = 3;

   int SpectralTrapezoidalAllL::truncationBwd(const int nN, const int l)
   {
      return std::max(nN - l/2, MIN_TRUNCATION);
   }

   int SpectralTrapezoidalAllL::index(const int i, const int k)
   {
      return i;
   }

   bool SpectralTrapezoidalAllL::isOptimal(const int nN, const int maxL)
   {
      return (this->truncationBwd(nN, maxL) > MIN_TRUNCATION);
   }

} // Tools
} // SpatialScheme
} // QuICC
