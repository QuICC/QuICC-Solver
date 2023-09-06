/** 
 * @file IMesher.cpp
 * @brief Source of the base for the spatial scheme mesher
 */

// System includes
//
#include <vector>

// Project includes
//
#include "QuICC/SpatialScheme/IMesher.hpp"

namespace QuICC {

namespace SpatialScheme {

   IMesher::IMesher(const GridPurpose::Id purpose)
      : mPurpose(purpose)
   {
   }

   void IMesher::init(const std::vector<int>& dims, const std::map<std::size_t,std::vector<std::size_t>>& options)
   {
      this->mDims = dims;
   }

} // SpatialScheme
} // QuICC
