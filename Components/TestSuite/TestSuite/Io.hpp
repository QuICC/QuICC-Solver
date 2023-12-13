/**
 * @file Io.hpp
 * @brief Useful I/O functions for tests
 */

#ifndef QUICC_TESTSUITE_IO_HPP
#define QUICC_TESTSUITE_IO_HPP

// System includes
//
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace TestSuite {

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


} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_IO_HPP
