/** 
 * @file ISphericalOperator.cpp
 * @brief Source of the base implementation of a spherical operator with y = ax + b and needing the harmonic degree
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/ISphericalOperator.hpp"

namespace QuICC {

namespace SparseSM {

namespace Chebyshev {

namespace LinearMap {

   ISphericalOperator::ISphericalOperator(const int rows, const int cols, const Scalar_t lower, const Scalar_t upper, const Scalar_t l)
      : ILinearMapOperator(rows, cols, lower, upper), mL(l)
   {
   }

   ISphericalOperator::Scalar_t ISphericalOperator::l() const
   {
      return this->mL;
   }

} // LinearMap
} // Chebyshev
} // SparseSM
} // QuICC
