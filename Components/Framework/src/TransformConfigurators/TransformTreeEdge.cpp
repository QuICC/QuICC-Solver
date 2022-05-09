/**
 * @file TransformTreeEdge.cpp
 * @brief Source of the implementation of the tranform tree edge
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/TransformConfigurators/TransformTree.hpp"

// Project includes
//
#include "QuICC/Arithmetics/None.hpp"
#include "QuICC/Arithmetics/Set.hpp"
#include "Profiler/Interface.hpp"

namespace QuICC {

namespace Transform {

   TransformTreeEdge::TransformTreeEdge(const std::size_t op, const int n)
      :mOpId(op), mN(n), mRecoverInput(false), mHoldInput(false), mRecoverOutId(-1), mCombinedOutId(-1), mFieldId(FieldType::SCALAR), mArithId(Arithmetics::Set::id()), mCombinedArithId(Arithmetics::None::id())
   {
      // Initialize output id
      this->mOutId.push_back(-1);
   }

   TransformTreeEdge::~TransformTreeEdge()
   {
   }

   std::size_t TransformTreeEdge::opId() const
   {
      return this->mOpId;
   }

   int TransformTreeEdge::nEdges(const int depth) const
   {
      if(depth == 0)
      {
         return this->mEdges.size();

      } else
      {
         int n = 0;
         for(auto it = this->mEdges.cbegin(); it != this->mEdges.cend(); ++it)
         {
            n += it->nEdges(depth - 1);
         }

         return n;
      }
   }

   TransformTreeEdge::EdgeType_crange TransformTreeEdge::edgeRange() const
   {
      Profiler::RegionFixture<4> fix("edgeRange");
      return std::make_pair(this->mEdges.begin(), this->mEdges.end());
   }

   TransformTreeEdge::EdgeType_range TransformTreeEdge::rEdgeRange()
   {
      return std::make_pair(this->mEdges.begin(), this->mEdges.end());
   }

   TransformTreeEdge& TransformTreeEdge::addEdge(const std::size_t op, const int n)
   {
      this->mEdges.push_back(TransformTreeEdge(op, n));

      return this->mEdges.back();
   }

   TransformTreeEdge::EdgeType_iterator TransformTreeEdge::delEdge(TransformTreeEdge::EdgeType_iterator edgeIt)
   {
      return this->mEdges.erase(edgeIt);
   }

   bool TransformTreeEdge::recoverInput() const
   {
      return this->mRecoverInput;
   }

   bool TransformTreeEdge::holdInput() const
   {
      return this->mHoldInput;
   }

   int TransformTreeEdge::recoverOutId() const
   {
      return this->mRecoverOutId;
   }

   int TransformTreeEdge::combinedOutId() const
   {
      return this->mCombinedOutId;
   }

   const std::vector<int>& TransformTreeEdge::outIds() const
   {
      return this->mOutId;
   }

   FieldType::Id TransformTreeEdge::fieldId() const
   {
      return this->mFieldId;
   }

   std::size_t TransformTreeEdge::arithId() const
   {
      return this->mArithId;
   }

   std::size_t TransformTreeEdge::combinedArithId() const
   {
      return this->mCombinedArithId;
   }

   void TransformTreeEdge::setInputInfo(const int recover, const int hold)
   {
      this->mRecoverInput = recover;

      this->mHoldInput = hold;
   }

   void TransformTreeEdge::setCombinedInfo(const std::size_t arithId, const int recoverId, const int holdId)
   {
      this->mCombinedArithId = arithId;

      this->mRecoverOutId = recoverId;

      this->mCombinedOutId = holdId;
   }

   void TransformTreeEdge::setArithId(const std::size_t arithId)
   {
      this->mArithId = arithId;
   }

   void TransformTreeEdge::setEnd(const int id, const FieldType::Id type, const std::size_t arith)
   {
      // Make sure output IDs is empty
      this->mOutId.clear();
      this->mOutId.push_back(id);

      this->mFieldId = type;

      this->mArithId = arith;
   }

   void TransformTreeEdge::setEnd(const std::vector<int>& id, const FieldType::Id type, const std::size_t arith)
   {
      // Make sure output IDs is empty
      this->mOutId.clear();
      this->mOutId = id;

      this->mFieldId = type;

      this->mArithId = arith;
   }

}
}
