/** 
 * @file IPlaneOperator.cpp
 * @brief Source of the base implementation of a plane layer operator with y = ax + b and k1, k2 Fourier modes
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/IPlaneOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   IPlaneOperator::IPlaneOperator(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t k1, const Scalar_t k2)
      : ILinearMapOperator(rows, cols, lower, upper), mK1(k1), mK2(k2)
   {
   }

   IPlaneOperator::Scalar_t IPlaneOperator::k1() const
   {
      return this->mK1;
   }

   IPlaneOperator::Scalar_t IPlaneOperator::k2() const
   {
      return this->mK2;
   }

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
