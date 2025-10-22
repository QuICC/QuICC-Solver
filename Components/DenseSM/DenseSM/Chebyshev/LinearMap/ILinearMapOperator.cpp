/**
 * @file ILinearMapOperator.cpp
 * @brief Source of the implementation of generic interface to a spherical shell
 * linear map dense operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ILinearMapOperator.hpp"
#include "include/QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ILinearMapOperator::ILinearMapOperator(const int rows, const int cols,
   const Scalar_t lower, const Scalar_t upper) :
    DenseSM::IMatrixSMOperator(rows, cols), mcLower(lower), mcUpper(upper)
{}

void ILinearMapOperator::computeQuadrature(Internal::Array& igrid,
   Internal::Array& iweights, const int size) const
{
   Polynomial::Quadrature::ChebyshevRule quad;
   quad.computeQuadrature(igrid, iweights, size, this->mcLower, this->mcUpper);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
