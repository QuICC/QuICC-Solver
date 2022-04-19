/**
 * @file Framework.cpp
 * @brief Implementation of Framework component version
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "QuICC/Version/Framework.hpp"

// Project includes
//

namespace QuICC {

namespace Version {

   const int Framework::MAJOR = QUICC_VERSION_FRAMEWORK_MAJOR;

   const int Framework::MINOR = QUICC_VERSION_FRAMEWORK_MINOR;

   const int Framework::PATCH = QUICC_VERSION_FRAMEWORK_PATCH;

   std::string Framework::version()
   {
      std::stringstream oss;

      oss << Framework::MAJOR << "." << Framework::MINOR << "." << Framework::PATCH;

      return oss.str();
   }

}
}
