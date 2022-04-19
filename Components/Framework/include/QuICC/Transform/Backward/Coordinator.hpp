/**
 * @file Coordinator.hpp
 * @brief Backward transform operator coordinator
 */

#ifndef QUICC_TRANSFORM_BACKWARD_COORDINATOR_HPP
#define QUICC_TRANSFORM_BACKWARD_COORDINATOR_HPP

// System includes
//
#include <string>
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Hasher.hpp"

namespace QuICC {

namespace Transform {

namespace Backward {

   /**
    * @brief Register new forward transform operator
    */
   template <typename TOperator> static std::size_t registerId(const std::string tag);

   template <typename TOperator> std::size_t registerId(const std::string s)
   {
      std::size_t id = Hasher::makeId(s);

      return id;
   }

}
}
}

#endif // QUICC_TRANSFORM_BACKWARD_COORDINATOR_HPP
