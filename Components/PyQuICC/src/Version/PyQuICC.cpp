/**
 * @file PyQuICC.cpp
 * @brief Implementation of PyQuICC component version
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "QuICC/Version/PyQuICC.hpp"

// Project includes
//

namespace QuICC {

namespace Version {

   const int PyQuICC::MAJOR = QUICC_VERSION_PYQUICC_MAJOR;

   const int PyQuICC::MINOR = QUICC_VERSION_PYQUICC_MINOR;

   const int PyQuICC::PATCH = QUICC_VERSION_PYQUICC_PATCH;

   std::string PyQuICC::version()
   {
      std::stringstream oss;

      oss << PyQuICC::MAJOR << "." << PyQuICC::MINOR << "." << PyQuICC::PATCH;

      return oss.str();
   }

}
}
