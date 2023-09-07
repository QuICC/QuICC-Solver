/** 
 * @file SpectralUniformAllL.cpp
 * @brief Source of the tools for uniform + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/SpectralUniformAllL.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int SpectralUniformAllL::truncationBwd(const int nN, const int j, const int k)
   {
      return nN;
   }

   int SpectralUniformAllL::index(const int i, const int j, const int k)
   {
      return i;
   }

   bool SpectralUniformAllL::isOptimal(const int nN, const int maxL)
   {
      return true;
   }

} // Tools
} // SpatialScheme
} // QuICC
