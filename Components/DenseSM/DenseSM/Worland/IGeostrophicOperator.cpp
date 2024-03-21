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
#include "Types/Internal/Math.hpp"
#include "IGeostrophicOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   IGeostrophicOperator::IGeostrophicOperator(const Scalar_t ugAlpha, const Scalar_t ugDBeta, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int q)
      : IEmbeddedOperator(rows, cols, alpha, dBeta), mcUgAlpha(ugAlpha), mcUgDBeta(ugDBeta)
   {
   }

   bool IGeostrophicOperator::isUgBasis() const
   {
      return this->isUgBasis(this->mcUgAlpha, this->mcUgDBeta);
   }

   bool IGeostrophicOperator::isUgBasis(const Scalar_t a, const Scalar_t b) const
   {
      bool isBasis = (a > -1 && b > -1);

      return isBasis;
   }

} // Worland
} // DenseSM
} // QuICC
