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
#include "QuICC/TransformConfigurators/TransformPath.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

   TransformPath::TransformPath(int startId, FieldType::Id fieldId)
      :mStartId(startId), mFieldId(fieldId)
   {
   }

   TransformPath::~TransformPath()
   {
   }

   void TransformPath::addEdge(const std::size_t opId, const int outId, const std::size_t arithId)
   {
      this->mEdges.push_back(TransformPathEdge(opId, outId, arithId));
   }

   void TransformPath::addEdge(const std::size_t opId, const std::pair<int,int>& outId, const std::size_t arithId)
   {
      this->mEdges.push_back(TransformPathEdge(opId, outId, arithId));
   }

   const TransformPathEdge& TransformPath::edge(const int i) const
   {
      return this->mEdges.at(i);
   }

   int TransformPath::nEdges() const
   {
      return this->mEdges.size();
   }

   int TransformPath::startId() const
   {
      return this->mStartId;
   }

   FieldType::Id TransformPath::fieldId() const
   {
      return this->mFieldId;
   }

}
}
