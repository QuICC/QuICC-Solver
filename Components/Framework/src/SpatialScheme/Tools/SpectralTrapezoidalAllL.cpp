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

   SpectralTrapezoidalAllL::SpectralTrapezoidalAllL()
      : IBaseAllL(static_cast<int>(MinimalTruncation::Triangular))
   {
   }

   SpectralTrapezoidalAllL::SpectralTrapezoidalAllL(const int min)
      : IBaseAllL(min)
   {
   }

   int SpectralTrapezoidalAllL::truncationBwd(const int nN, const int j, const int k)
   {
      return std::max(nN - j/2, this->min());
   }

   int SpectralTrapezoidalAllL::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool SpectralTrapezoidalAllL::isOptimal(const int nN, const int maxL)
   {
      return (nN - maxL/2 >= this->min());
   }

} // Tools
} // SpatialScheme
} // QuICC
