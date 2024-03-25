/**
 * @file IGeostrophicOperator.cpp
 * @brief Source of the implementation of the base for a geostrophic projection operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "IGeostrophicOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   IGeostrophicOperator::IGeostrophicOperator(const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const int rows, const int cols)
      : IWorlandOperator(rows, cols), mcUgAlpha(ugAlpha), mcUgBeta(ugBeta), mcIsGenericBasis(isGenericBasis)
   {
   }

   IGeostrophicOperator::IGeostrophicOperator(const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta)
      : IWorlandOperator(rows, cols, alpha, dBeta), mcUgAlpha(ugAlpha), mcUgBeta(ugBeta), mcIsGenericBasis(isGenericBasis)
   {
   }

} // Worland
} // DenseSM
} // QuICC
