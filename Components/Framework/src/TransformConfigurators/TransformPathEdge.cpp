/**
 * @file TransformTree.cpp
 * @brief Source of the implementation of the tranform tree
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/TransformConfigurators/TransformPathEdge.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   TransformPathEdge::TransformPathEdge(const std::size_t opId, const int outId, const std::size_t arithId)
      :mOpId(opId), mArithId(arithId)
   {
      this->mOutId.push_back(outId);
   }

   TransformPathEdge::TransformPathEdge(const std::size_t opId, const std::pair<int,int>& outId, const std::size_t arithId)
      :mOpId(opId), mArithId(arithId)
   {
      this->mOutId.push_back(outId.first);
      this->mOutId.push_back(outId.second);
   }

   TransformPathEdge::~TransformPathEdge()
   {
   }

   std::size_t TransformPathEdge::opId() const
   {
      return this->mOpId;
   }

   const std::vector<int>& TransformPathEdge::outId() const
   {
      return this->mOutId;
   }

   std::size_t TransformPathEdge::arithId() const
   {
      return this->mArithId;
   }

}
}
