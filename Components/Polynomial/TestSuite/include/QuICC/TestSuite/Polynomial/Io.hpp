/**
 * @file Io.hpp
 * @brief Generic IO functions for tests
 */

#ifndef QUICC_TESTSUITE_POLYNOMIAL_IO_HPP
#define QUICC_TESTSUITE_POLYNOMIAL_IO_HPP

// Configuration includes
//

// System includes
//
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdexcept>

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace TestSuite {

namespace Polynomial {

   /**
    * @brief Write real data to file
    */
   void writeData(const std::string& path, const Matrix& outData);

   /**
    * @brief Write complex data to file
    */
   void writeData(const std::string& path, const MatrixZ& outData);

   /**
    * @brief Read real data from file
    */
   void readData(Matrix& inData, const std::string& path);

   /**
    * @brief Read complex data from file
    */
   void readData(MatrixZ& inData, const std::string& path);

}
}
}

#endif //QUICC_TESTSUITE_POLYNOMIAL_IO_HPP
