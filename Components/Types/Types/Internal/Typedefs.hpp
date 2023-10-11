/**
 * @file Typedefs.hpp
 * @brief Definition of typedefs for internal computations
 */

#ifndef QUICC_TYPES_INTERNAL_TYPEDEFS_HPP
#define QUICC_TYPES_INTERNAL_TYPEDEFS_HPP

// System includes
//
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"

namespace QuICC {
namespace Internal {

/// @brief Typedef for an array of internal float values
typedef Eigen::Matrix<MHDFloat, Eigen::Dynamic, 1> Array;
/// @brief Typedef for a coefficient wise array of float values
typedef Eigen::Array<MHDFloat, Eigen::Dynamic, 1> ACoeff;
/// @brief Typedef for a matrix of internal float values
typedef Eigen::Matrix<MHDFloat, Eigen::Dynamic, Eigen::Dynamic> Matrix;
/// @brief Typedef for a sparse matrix of internal float values
typedef Eigen::SparseMatrix<MHDFloat> SparseMatrix;
/// @brief Typedef for the internal Array type of extended float values
typedef Eigen::Matrix<MHDLong, Eigen::Dynamic, 1> ArrayL;
/// @brief Typedef for the internal Matrix type of extended float values
typedef Eigen::Matrix<MHDLong, Eigen::Dynamic, Eigen::Dynamic> MatrixL;
/// @brief Typedef for an smart reference counting pointer on an array of
/// internal real values
typedef std::shared_ptr<Array> SharedArray;
/// @brief Typedef for an smart reference counting pointer on an matrix of
/// internal real values
typedef std::shared_ptr<Matrix> SharedMatrix;

} // namespace Internal
} // namespace QuICC

#endif // QUICC_TYPES_INTERNAL_TYPEDEFS_HPP
