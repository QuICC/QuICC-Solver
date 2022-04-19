/** 
 * @file ICosts.cpp
 * @brief Source of the base for a cost based scheme implementations
 */

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/SpatialScheme/ICosts.hpp"

// Project includes
//

namespace QuICC {

namespace SpatialScheme {

   ICosts::ICosts(const int dims)
      : mCosts(dims), mScalings(dims), mMemory(dims)
   {
   }

   ICosts::~ICosts()
   {
   }

   void ICosts::resetLoad()
   {
      // Clear load list
      this->mLoadList.clear();

      // Clear optimal loads
      while (!this->mOptimalLoad.empty())
      {
         this->mOptimalLoad.pop();
      }

      // Clear load sums
      this->mLoad.clear();
      
      // Clear regularized loads per part
      this->mRegularLoad.clear();
   }

   void ICosts::updateLoad(const int parts)
   {
      // Create typedef to simplify notation
      typedef  std::multimap<int, int>::const_iterator MapIt;

      // Const iterator to map object
      MapIt it;

      // Pair of iterator to work as range
      std::pair<MapIt, MapIt> range;

      // Loop over all parts
      for(int i = 0; i < parts; ++i)
      {
         // Reset loads
         this->mLoad.at(i) = 0;

         // Get assigned loads
         range = this->mRegularLoad.equal_range(i);

         // Loop over the assigned loads
         for(it = range.first; it != range.second; ++it)
         {
            this->mLoad.at(i) += it->second;
         }
      }
   }

   Array ICosts::loadWeights()
   {
      return this->mCosts.array() * this->mScalings.array();
   }

   double ICosts::memoryScore(SharedResolution spRes)
   {
      #ifdef QUICC_MEMORYUSAGE_HIGH
         return 1.0;
      #else
         return this->mMemory.prod();
      #endif //QUICC_MEMORYUSAGE_HIGH
   }

   void ICosts::setCost(MHDFloat c, Dimensions::Transform::Id id)
   {
      // Assert for positive cost
      assert(c > 0);

      // Assert for number of transforms
      assert(static_cast<int>(id) < this->mCosts.size());

      this->mCosts(static_cast<int>(id)) = c;
   }

   void ICosts::setScaling(MHDFloat c, Dimensions::Transform::Id id)
   {
      // Assert for positive cost
      assert(c > 0);

      // Assert for number of transforms
      assert(static_cast<int>(id) < this->mScalings.size());

      this->mScalings(static_cast<int>(id)) = c;
   }

   void ICosts::setMemory(MHDFloat c, Dimensions::Transform::Id id)
   {
      // Assert for positive cost
      assert(c > 0);

      // Assert for number of transforms
      assert(static_cast<int>(id) < this->mMemory.size());

      this->mMemory(static_cast<int>(id)) = c;
   }

}
}
