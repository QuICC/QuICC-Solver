/**
 * @file RadialPowerDivY1D1Y1.cpp
 * @brief Source of the implementation of the Chebyshev radial power 1/Y^1 D^1 Y^1
 * reductor, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/RadialPowerDivY1D1Y1.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

void RadialPowerDivY1D1Y1::initOperator() const
{
      // Check for division by 0!
      assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

      ::QuICC::SparseSM::Chebyshev::LinearMap::I1 op(this->mspSetup->specSize()+2,this->mspSetup->specSize()+2, this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.solver().setOperator(op.mat(), 1, 1);

      ::QuICC::SparseSM::Chebyshev::LinearMap::Y1 opY1(this->mspSetup->specSize()+1,this->mspSetup->specSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.solver().setSpectralOperator(opY1.mat(), 1);

      Internal::Array igrid, iweights;
      Polynomial::Quadrature::ChebyshevRule quad;
      quad.computeQuadrature(igrid, iweights, 2*this->mspSetup->fwdSize(), this->mspSetup->lower(), this->mspSetup->upper());
      this->mBackend.setScaler(igrid.array().pow(-1).cast<MHDFloat>().matrix());
}

void RadialPowerDivY1D1Y1::initBackend() const
{
   // Call parent initializer
   ILinearMapRadialPower::initBackend();

   // Initialize the solver
   this->mBackend.addSolver(1);
}

void RadialPowerDivY1D1Y1::applyPreOperator(Matrix& tmp, const Matrix& in) const
{
   this->mBackend.input(tmp, in);
   auto specOp = this->mBackend.solver().getSpectralOperator();
   tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
   this->mBackend.getSolution(tmp, 1, 1);
}

void RadialPowerDivY1D1Y1::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
{
   assert(rOut.cols() == 1);
   this->mBackend.outputGrid(rOut, tmp);
}

void RadialPowerDivY1D1Y1::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
   auto specOp = this->mBackend.solver().getSpectralOperator();
   tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
   this->mBackend.getSolution(tmp, 1, 1);
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
