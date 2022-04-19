/** 
 * @file IOperator.cpp
 * @brief Source of forward projection operator IOperator 
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Forward/IOperator.hpp"

// Project includes
//
#include "QuICC/Transform/Forward/Coordinator.hpp"

namespace QuICC {

namespace Transform {

namespace Forward {

   IOperator::IOperator(const std::string tag)
      : mTag(tag)
   {
   }

   IOperator::~IOperator()
   {
   }

   std::string IOperator::tag() const
   {
      return this->mTag;
   }

}
}
}
