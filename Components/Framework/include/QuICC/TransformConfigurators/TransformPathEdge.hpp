/**
 * @file TransformPathEdge.hpp
 * @brief This class defines an edge in transform path
 */

#ifndef QUICC_TRANSFORM_TRANSFORMPATHEDGE_HPP
#define QUICC_TRANSFORM_TRANSFORMPATHEDGE_HPP

// Configuration includes
//

// System includes
//
#include <set>
#include <vector>

// External includes
//

// Project includes
//

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines an edge in a transform paht
    */
   class TransformPathEdge
   {
      public:
         /**
          * @brief Contructor for edge with single ID for output field
          */
         TransformPathEdge(const std::size_t opId, const int outId, const std::size_t arithId);

         /**
          * @brief Contructor for edge with ID pair for output field
          */
         TransformPathEdge(const std::size_t opId, const std::pair<int,int>& outId, const std::size_t arithId);

         /**
          * @brief Destructor
          */
         ~TransformPathEdge();

         /**
          * @brief Get ID of operator
          */
         std::size_t opId() const;

         /**
          * @brief Get IDs of output field
          */
         const std::vector<int>& outId() const;

         /**
          * @brief Get ID of arithmetic operation
          */
         std::size_t arithId() const;

      private:
         /**
          * @brief ID of operator
          */
         std::size_t mOpId;

         /**
          * @brief ID of the arithmetic operation
          */
         std::size_t mArithId;

         /**
          * @brief IDs of the output field
          */
         std::vector<int>  mOutId;
   };

}
}

#endif // QUICC_TRANSFORM_TRANSFORMPATHEDGE_HPP
