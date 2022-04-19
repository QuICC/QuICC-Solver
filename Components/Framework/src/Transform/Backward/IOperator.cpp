/** 
 * @file IOperator.cpp
 * @brief Source of the generic backward transform operator
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Transform/Backward/IOperator.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace Backward {

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
