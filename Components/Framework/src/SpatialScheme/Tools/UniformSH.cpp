/** 
 * @file UniformSH.cpp
 * @brief Source of the tools for uniform + SH spatial schemes
 */

// System includes
//
#include <set>

// Project includes
//
#include "QuICC/SpatialScheme/Tools/UniformSH.hpp"

namespace QuICC {

namespace SpatialScheme {

namespace Tools {

   int UniformSH::truncationFwd(const int nN, const int l)
   {
      return nN;
   }

   int UniformSH::truncationBwd(const int nN, const int l)
   {
      return nN;
   }

   int UniformSH::index(const int i, const int k)
   {
      return i;
   }

   bool UniformSH::isOptimal(const int nN, const int maxL)
   {
      return true;
   }

} // Tools
} // SpatialScheme
} // QuICC
