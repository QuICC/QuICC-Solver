/**
 * @file Io.hpp
 * @brief Useful I/O functions for tests
 */

#ifndef QUICC_TESTSUITE_TRANSFORM_IO_HPP
#define QUICC_TESTSUITE_TRANSFORM_IO_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {

namespace TestSuite {

namespace Transform {

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
   void readList(Array& inData, const std::string& path);

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

#endif // QUICC_TESTSUITE_TRANSFORM_IO_HPP
