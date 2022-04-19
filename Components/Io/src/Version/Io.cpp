/**
 * @file Io.cpp
 * @brief Implementation of Io component version
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "QuICC/Version/Io.hpp"

// Project includes
//

namespace QuICC {

namespace Version {

   const int Io::MAJOR = QUICC_VERSION_IO_MAJOR;

   const int Io::MINOR = QUICC_VERSION_IO_MINOR;

   const int Io::PATCH = QUICC_VERSION_IO_PATCH;

   std::string Io::version()
   {
      std::stringstream oss;

      oss << Io::MAJOR << "." << Io::MINOR << "." << Io::PATCH;

      return oss.str();
   }

}
}
