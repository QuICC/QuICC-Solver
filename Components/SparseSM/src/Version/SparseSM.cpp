/**
 * @file SparseSM.cpp
 * @brief Implementation of SparseSM component version
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "QuICC/Version/SparseSM.hpp"

// Project includes
//

namespace QuICC {

namespace Version {

   const int SparseSM::MAJOR = QUICC_VERSION_SPARSESM_MAJOR;

   const int SparseSM::MINOR = QUICC_VERSION_SPARSESM_MINOR;

   const int SparseSM::PATCH = QUICC_VERSION_SPARSESM_PATCH;

   std::string SparseSM::version()
   {
      std::stringstream oss;

      oss << SparseSM::MAJOR << "." << SparseSM::MINOR << "." << SparseSM::PATCH;

      return oss.str();
   }

}
}
