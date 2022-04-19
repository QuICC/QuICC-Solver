/**
 * @file TransformPath.hpp
 * @brief This class defines a general transform path
 */

#ifndef QUICC_TRANSFORM_TRANSFORMPATH_HPP
#define QUICC_TRANSFORM_TRANSFORMPATH_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Arithmetics/Set.hpp"
#include "QuICC/TransformConfigurators/TransformPathEdge.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class describes a backward transform tree branch
    */
   class TransformPath
   {
      public:
         /**
          * @brief Contructor for transform path
          *
          * @param startId Field component ID
          * @param fieldId Field type ID
          */
         TransformPath(int startId, FieldType::Id fieldId);

         /**
          * @brief Destructor
          */
         ~TransformPath();

         /**
          * @brief Add an edge to the tranform path with single ID for output field
          */
         void addEdge(const std::size_t opId, const int outId = -1, const std::size_t arithId = Arithmetics::Set::id());

         /**
          * @brief Add an edge to the transform path with ID pair for output field
          */
         void addEdge(const std::size_t opId, const std::pair<int,int>& outId, const std::size_t arithId);

         /**
          * @brief Get edge
          */
         const TransformPathEdge& edge(const int i) const;

         /**
          * @brief Get the number of edges
          */
         int nEdges() const;

         /**
          * @brief Get starting component ID
          */
         int startId() const;

         /**
          * @brief Get field type ID
          */
         FieldType::Id fieldId() const;

      private:
         /**
          * @brief Starting component for transform path
          */
         int mStartId;

         /**
          * @brief Field type required for transform branch
          */
         FieldType::Id mFieldId;

         /**
          * @brief Vector of edges of the branch
          */
         std::vector<TransformPathEdge> mEdges;
   };

}
}

#endif // QUICC_TRANSFORM_TRANSFORMPATH_HPP
