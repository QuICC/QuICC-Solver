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
/// @brief namespace for TestSuite common utilities
namespace TestSuite {

/// @brief Write real data to file
/// @param path
/// @param outData dense real matrix
void writeData(const std::string& path, const Matrix& outData);

/// @brief Write complex data to file
/// @param path
/// @param outData dense complex matrix
void writeData(const std::string& path, const MatrixZ& outData);

/// @brief Write real data to file
/// @param path
/// @param outData sparse real matrix
void writeData(const std::string& path, const SparseMatrix& outData);

/// @brief
/// @param inData
/// @param path
void readList(Array& inData, const std::string& path);

/// @brief Read real data from file
/// @param inData dense real matrix
/// @param path
void readData(Matrix& inData, const std::string& path);

/// @brief Read complex data from file
/// @param inData dense complex matrix
/// @param path
void readData(MatrixZ& inData, const std::string& path);

/// @brief Read real data from file
/// @param inData sparse real matrix
/// @param path
void readData(SparseMatrix& inData, const std::string& path);


} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_IO_HPP
