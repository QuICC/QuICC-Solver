/**
 * @file TransformTreeEdge.hpp
 * @brief This class defines a general transform tree edge
 */

#ifndef QUICC_TRANSFORM_TRANSFORMTREEEDGE_HPP
#define QUICC_TRANSFORM_TRANSFORMTREEEDGE_HPP

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"

namespace QuICC {

namespace Transform {

   // Forward declaration
   class TransformTreeEdge;

   /**
    * @brief Specialisation to end recursion
    */
   class TransformTreeEdge
   {
      public:
         /// Useful typedefs
         typedef std::vector<TransformTreeEdge>::const_iterator EdgeType_citerator;
         typedef std::pair<EdgeType_citerator,EdgeType_citerator> EdgeType_crange;
         typedef std::vector<TransformTreeEdge>::iterator EdgeType_iterator;
         typedef std::pair<EdgeType_iterator,EdgeType_iterator> EdgeType_range;

         /**
          * @brief Contructor for operation
          */
         TransformTreeEdge(const std::size_t op, const int n);

         /**
          * @brief Destructor
          */
         ~TransformTreeEdge();

         /**
          * @brief Get operator ID
          */
         std::size_t opId() const;

         /**
          * @brief Get the number of edges
          */
         int nEdges(const int depth) const;

         /**
          * @brief Get vector ranges of edges
          */
         EdgeType_crange edgeRange() const;

         /**
          * @brief Get vector ranges of edges
          */
         EdgeType_range rEdgeRange();

         /**
          * @brief Add a edge to tree
          */
         TransformTreeEdge& addEdge(const std::size_t op, const int n);

         /**
          * @brief Delete an edge to tree
          */
         EdgeType_iterator delEdge(EdgeType_iterator edgeIt);

         /**
          * @brief Get out component id
          */
         template <typename TId> TId outId(const int i = 0) const;

         /**
          * @brief Get all output component ids
          */
         const std::vector<int>& outIds() const;

         /**
          * @brief Recover input from other calculations
          */
         bool recoverInput() const;

         /**
          * @brief Hold input for other calculations
          */
         bool holdInput() const;

         /**
          * @brief Recover output from other calculations
          */
         int recoverOutId() const;

         /**
          * @brief Output ID to hold combined calculation
          */
         int combinedOutId() const;

         /**
          * @brief Get the field type
          */
         FieldType::Id fieldId() const;

         /**
          * @brief Get the arithmetic operation
          */
         std::size_t arithId() const;

         /**
          * @brief Get the arithmetic operation on the combined field
          */
         std::size_t combinedArithId() const;

         /**
          * @brief Set the arithmetic operation
          */
         void setArithId(const std::size_t);

         /**
          * @brief Set recovery and hold information for input data
          */
         void setInputInfo(const int recover, const int hold);

         /**
          * @brief Set output recovery and combination inforation
          */
         void setCombinedInfo(const std::size_t, const int recoverId, const int holdId);

         /**
          * @brief Set end component as vector of IDs
          */
         void setEnd(const std::vector<int>& id, const FieldType::Id type, const std::size_t arith);

         /**
          * @brief Set end component as a single ID
          */
         void setEnd(const int id, const FieldType::Id type, const std::size_t arith);

      private:
         /**
          * @brief Operation attached to edge
          */
         std::size_t  mOpId;

         /**
          * @brief Multiplicity of edge (used by how many branches)
          */
         int mN;

         /**
          * @brief Recover input from other computations
          */
         bool mRecoverInput;

         /**
          * @brief Hold input for other computations
          */
         bool mHoldInput;

         /**
          * @brief Recover output from other computations (ID >= 0 if recovery is necessary)
          */
         int mRecoverOutId;

         /**
          * @brief Hold combined output for other computations (ID >= 0 if hold is necessary)
          */
         int mCombinedOutId;

         /**
          * @brief Field type
          */
         FieldType::Id  mFieldId;

         /**
          * @brief Arithmetic operation to store result
          */
         std::size_t   mArithId;

         /**
          * @brief Arithmetic operation for combined output
          */
         std::size_t   mCombinedArithId;

         /**
          * @brief Output Id components
          */
         std::vector<int> mOutId;

         /**
          * @brief Vector of connected egdes
          */
         std::vector<TransformTreeEdge>  mEdges;
   };

   template <typename TId> inline TId TransformTreeEdge::outId(const int i) const
   {
      return static_cast<TId>(this->mOutId.at(i));
   }
}
}

#endif // QUICC_TRANSFORM_TRANSFORMTREEEDGE_HPP
