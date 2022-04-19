/** 
 * @file CoreResolution.cpp
 * @brief Source of the resolution object for a single CPU
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "QuICC/Resolutions/CoreResolution.hpp"

// Project includes
//

namespace QuICC {

   CoreResolution::CoreResolution(const std::vector<SharedTransformResolution>& transformRes)
      : mTransforms(transformRes)
   {
   }

   CoreResolution::~CoreResolution()
   {
   }

   int CoreResolution::nDim() const
   {
      return this->mTransforms.size();
   }

   bool CoreResolution::isCleared() const
   {
      bool cleared = true;
      for(const auto& r: this->mTransforms)
      {
         cleared = cleared && r->isCleared();
      }
      return cleared;
   }

   SharedCTransformResolution CoreResolution::dim(const Dimensions::Transform::Id id) const
   {
      // Check for correct sizes
      assert(static_cast<int>(id) >= 0);
      assert(static_cast<size_t>(id) < this->mTransforms.size());

      return this->mTransforms.at(static_cast<int>(id));
   }

}
