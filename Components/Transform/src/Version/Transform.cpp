/**
 * @file Transform.cpp
 * @brief Implementation of Transform component version
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "QuICC/Version/Transform.hpp"

// Project includes
//

namespace QuICC {

namespace Version {

   const int Transform::MAJOR = QUICC_VERSION_TRANSFORM_MAJOR;

   const int Transform::MINOR = QUICC_VERSION_TRANSFORM_MINOR;

   const int Transform::PATCH = QUICC_VERSION_TRANSFORM_PATCH;

   std::string Transform::version()
   {
      std::stringstream oss;

      oss << Transform::MAJOR << "." << Transform::MINOR << "." << Transform::PATCH;

      return oss.str();
   }

}
}
