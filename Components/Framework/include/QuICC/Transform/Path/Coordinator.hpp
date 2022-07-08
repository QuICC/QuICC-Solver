/**
 * @file Coordinator.hpp
 * @brief Transform path IDcoordinator
 */

#ifndef QUICC_TRANSFORM_PATH_COORDINATOR_HPP
#define QUICC_TRANSFORM_PATH_COORDINATOR_HPP

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

namespace Path {

   /**
    * @brief Register new ID
    */
   template <typename TId> static std::size_t registerId(const std::string tag);

   template <typename TId> std::size_t registerId(const std::string s)
   {
      std::size_t id = Hasher::makeId(s);

      return id;
   }

}
}
}

#endif // QUICC_TRANSFORM_PATH_COORDINATOR_HPP
