/**
 * @file Typedefs.hpp
 * @brief Definition of typedefs for muliple precision computations
 */

#ifndef QUICC_TYPES_MP_TYPEDEFS_HPP
#define QUICC_TYPES_MP_TYPEDEFS_HPP

// System includes
//
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

// Project includes
//
#include "Types/MP/BasicTypes.hpp"

namespace QuICC {
namespace MP {

/**
 * @name Array types typedefs
 */
//@{
/// Typedef for an array of multiple precision float values
typedef Eigen::Matrix<MHDFloat, Eigen::Dynamic, 1> Array;
//@}

/**
 * @name Coefficient wise array types typedefs
 */
//@{
/// Typedef for an array of float values
typedef Eigen::Array<MHDFloat, Eigen::Dynamic, 1> ACoeff;
//@}

/**
 * @name Matrix types typedefs
 */
//@{
/// Typedef for a matrix of multiple precision float values
typedef Eigen::Matrix<MHDFloat, Eigen::Dynamic, Eigen::Dynamic> Matrix;
//@}
//
/**
 * @name Sparse Matrix types typedefs
 */
//@{
/// Typedef for a sparse matrix of float values
typedef Eigen::SparseMatrix<MHDFloat> SparseMatrix;
//@}

/**
 * @name Shared pointer typedefs
 */
//@{
/// Typedef for an smart reference counting pointer on an array of multiple
/// precision real values
typedef std::shared_ptr<Array> SharedArray;
/// Typedef for an smart reference counting pointer on an matrix of real values
typedef std::shared_ptr<Matrix> SharedMatrix;
//@}

} // namespace MP
} // namespace QuICC

#endif // QUICC_TYPES_MP_TYPEDEFS_HPP
