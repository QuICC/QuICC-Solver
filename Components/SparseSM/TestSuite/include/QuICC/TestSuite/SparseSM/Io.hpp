/**
 * @file Io.hpp
 * @brief Generic IO functions for tests
 */

#ifndef QUICC_TESTSUITE_SPARSESM_IO_HPP
#define QUICC_TESTSUITE_SPARSESM_IO_HPP

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
#include "QuICC/Typedefs.hpp"
#include "QuICC/Precision.hpp"

namespace QuICC {

namespace TestSuite {

namespace SparseSM {

   /**
    * @brief Write real data to file
    */
   void writeData(const std::string& path, const SparseMatrix& outData);

   /**
    * @brief Write real data to file
    */
   void writeData(const std::string& path, const Matrix& outData);

   /**
    * @brief Read real data from file
    */
   void readData(SparseMatrix& inData, const std::string& path);

   /**
    * @brief Read real data from file
    */
   void readData(Matrix& inData, const std::string& path);

   /**
    * @brief Read real data from file
    */
   void readList(Array& inData, const std::string& path);

}
}
}

#endif //QUICC_TESTSUITE_SPARSESM_IO_HPP
