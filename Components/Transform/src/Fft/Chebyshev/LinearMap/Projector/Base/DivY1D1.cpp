/**
 * @file DivY1D1.cpp
 * @brief Source of the implementation of the Chebyshev 1/Y D projector, with linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/DivY1D1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

void DivY1D1<base_t>::initOperator() const
{
   // Check for division by 0!
      assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);
      // set up operator for f = I1 g as in D.hpp
      ::QuICC::SparseSM::Chebyshev::LinearMap::I1 op(this->mspSetup->specSize()+1,
                                                     this->mspSetup->specSize()+1,
                                                     this->mspSetup->lower(), 
                                                     this->mspSetup->upper());
      this->mBackend.solver().setOperator(op.mat());

      // set up scaler for division by y
      Internal::Array igrid, iweights;
      Polynomial::Quadrature::ChebyshevRule quad;
      quad.computeQuadrature(igrid, iweights, this->mspSetup->fwdSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setScaler(igrid.array().pow(-1).cast<MHDFloat>().matrix());
}

void DivY1D1<base_t>::initBackend() const
{
   // Call parent initializer
   ILinearMapProjector::initBackend();

   // Initialize the solver
   this->mBackend.addSolver();
}

void DivY1D1<base_t>::applyPreOperator(Matrix& tmp, const Matrix& in) const
{
   this->mBackend.input(tmp, in, 1);

   this->mBackend.getSolution(tmp, 1);
}

void DivY1D1<base_t>::applyPostOperator(Matrix& rOut) const 
{
   this->mBackend.outputScale(rOut);
}

void DivY1D1<base_t>::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, 1, useReal);
   this->mBackend.getSolution(tmp, 1);
}

void DivY1D1<base_t>::applyPostOperator(MatrixZ& rOut, const Matrix& tmp,
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
