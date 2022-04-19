/**
 * @file Polynomial.cpp
 * @brief Implementation of Polynomial component version
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "QuICC/Version/Polynomial.hpp"

// Project includes
//

namespace QuICC {

namespace Version {

   const int Polynomial::MAJOR = QUICC_VERSION_POLYNOMIAL_MAJOR;

   const int Polynomial::MINOR = QUICC_VERSION_POLYNOMIAL_MINOR;

   const int Polynomial::PATCH = QUICC_VERSION_POLYNOMIAL_PATCH;

   std::string Polynomial::version()
   {
      std::stringstream oss;

      oss << Polynomial::MAJOR << "." << Polynomial::MINOR << "." << Polynomial::PATCH;

      return oss.str();
   }

}
}
