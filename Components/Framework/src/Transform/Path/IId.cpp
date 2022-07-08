/** 
 * @file IId.cpp
 * @brief Source of transform path ID
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Path/IId.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Path {

   IId::IId(const std::string tag)
      : mTag(tag)
   {
   }

   IId::~IId()
   {
   }

   std::string IId::tag() const
   {
      return this->mTag;
   }

}
}
}
