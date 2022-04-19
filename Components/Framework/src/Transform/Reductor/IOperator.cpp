/** 
 * @file IOperator.cpp
 * @brief Source of the generic reductor transform operator
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Reductor/IOperator.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Reductor {

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
