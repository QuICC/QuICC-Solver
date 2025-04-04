/**
 * @file DivY2.cpp
 * @brief Source of the implementation of the Chebyshev 1/R^2 projector
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/DivY2.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

void DivY2<base_t>::initOperator() const
{
   // Check for division by 0!
   assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

   Internal::Array igrid, iweights;
   Polynomial::Quadrature::ChebyshevRule quad;
   quad.computeQuadrature(igrid, iweights, this->mspSetup->fwdSize(),
      this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.setScaler(igrid.array().pow(-2).cast<MHDFloat>().matrix());
}

void DivY2<base_t>::applyPreOperator(Matrix& tmp, const Matrix& in) const
{
   this->mBackend.input(tmp, in);
}

void DivY2<base_t>::applyPostOperator(Matrix& rOut) const
{
   this->mBackend.outputScale(rOut);
}

void DivY2<base_t>::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
}

void DivY2<base_t>::applyPostOperator(MatrixZ& rOut, const Matrix& tmp,
   const bool useReal) const
{
   this->mBackend.outputScale(rOut, tmp, useReal);
}

} // namespace Projector
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
