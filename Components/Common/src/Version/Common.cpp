/**
 * @file Common.cpp
 * @brief Implementation of Common component version
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "QuICC/Version/Common.hpp"

// Project includes
//

namespace QuICC {

namespace Version {

   const int Common::MAJOR = QUICC_VERSION_COMMON_MAJOR;

   const int Common::MINOR = QUICC_VERSION_COMMON_MINOR;

   const int Common::PATCH = QUICC_VERSION_COMMON_PATCH;

   std::string Common::version()
   {
      std::stringstream oss;

      oss << Common::MAJOR << "." << Common::MINOR << "." << Common::PATCH;

      return oss.str();
   }

}
}
