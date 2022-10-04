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
}
}

#endif // QUICC_DIMENSIONS_DIMENSIONTOOLS_HPP
