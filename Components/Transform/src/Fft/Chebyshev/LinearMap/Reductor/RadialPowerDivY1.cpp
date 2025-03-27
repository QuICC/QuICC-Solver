/**
 * @file RadialPowerDivY1.cpp
 * @brief Source of the implementation of the Chebyshev radial power 1/Y^1
 * reductor, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/RadialPowerDivY1.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

void RadialPowerDivY1::initOperator() const
{
   // Check for division by 0!
   assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

   Internal::Array igrid, iweights;
   Polynomial::Quadrature::ChebyshevRule quad;
   quad.computeQuadrature(igrid, iweights, 2*this->mspSetup->fwdSize(), this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.setScaler(igrid.array().pow(-1).cast<MHDFloat>().matrix());
}

void RadialPowerDivY1::applyPreOperator(Matrix& tmp, const Matrix& in) const
{
   this->mBackend.input(tmp, in);
}

void RadialPowerDivY1::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
{
   assert(rOut.cols() == 1);
   this->mBackend.outputGrid(rOut, tmp);
}

void RadialPowerDivY1::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
