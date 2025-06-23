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
/// @param path path to the file
/// @param outData dense real matrix
void writeData(const std::string& path, const Matrix& outData);

/// @brief Write complex data to file
/// @param path path to the file
/// @param outData dense complex matrix
void writeData(const std::string& path, const MatrixZ& outData);

/// @brief Write real data to file
/// @param path path to the file
/// @param outData sparse real matrix
void writeData(const std::string& path, const SparseMatrix& outData);

/// @brief Read list of values (Eigen array)
/// @param inData storage for values
/// @param path path to the file
void readList(Array& inData, const std::string& path);

/// @brief Read list of values (std::vector)
/// @param inData storage for values
/// @param path path to the file
void readList(std::vector<MHDFloat>& inData, const std::string& path);

/// @brief Read real data from file
/// @param inData dense real matrix
/// @param path path to the file
void readData(Matrix& inData, const std::string& path);

/// @brief Read real data separated in blocks from file
/// @param inData dense real matrix
/// @param path path to the file
void readBlockData(std::vector<Matrix>& inData, const std::string& path);

/// @brief Read complex data from file
/// @param inData dense complex matrix
/// @param path path to the file
void readData(MatrixZ& inData, const std::string& path);

/// @brief Read real data from file
/// @param inData sparse real matrix
/// @param path path to the file
void readData(SparseMatrix& inData, const std::string& path);

/// @brief Read lines from file
/// @param lines storage for lines
/// @param path path to the file
void readLines(std::vector<std::string>& lines, const std::string& path);

/// @brief Split line into words
/// @param words storage for words
/// @param line	line to split
/// @param delim delimiter to use for splitting
void splitLine(std::vector<std::string>& words, const std::string& line, const char delim);

void getCommand(int& argc, std::vector<char *> argv, const std::vector<std::string>& options);

} // namespace TestSuite
} // namespace QuICC

#endif // QUICC_TESTSUITE_IO_HPP
