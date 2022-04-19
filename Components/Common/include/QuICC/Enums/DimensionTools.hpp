/**
 * @file DimensionTools.hpp
 * @brief Definition of some useful tools for the dimensions enums 
 */

#ifndef QUICC_DIMENSIONS_DIMENSIONTOOLS_HPP
#define QUICC_DIMENSIONS_DIMENSIONTOOLS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"

namespace QuICC {

namespace Dimensions {

      /**
       * @brief Jump to another transform dimension at runtime
       *
       * @param id      Base id
       * @param step    Size of the dimension jump
       */
      Transform::Id  jump(Transform::Id id, int step);

      /**
       * @brief Jump to another dimension at compilation time
       */
      template<Transform::Id TID, int STEP> struct TransformJump
      {
         constexpr static const Transform::Id id = static_cast<Transform::Id>(static_cast<int>(TID)+STEP);
      };
}
}

#endif // QUICC_DIMENSIONS_DIMENSIONTOOLS_HPP
