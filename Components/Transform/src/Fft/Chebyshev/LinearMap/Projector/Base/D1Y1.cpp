/**
 * @file D1Y1.cpp
 * @brief Source of the implementation of the Chebyshev D Y projector, with
 * linear map y = ax + b
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/Base/D1Y1.hpp"

namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Projector {

void D1Y1<base_t>::initOperator() const
{
   ::QuICC::SparseSM::Chebyshev::LinearMap::I1 op(
      this->mspSetup->specSize() + 2, this->mspSetup->specSize() + 2,
      this->mspSetup->lower(), this->mspSetup->upper());

   this->mBackend.solver().setOperator(op.mat(), 1, 1);

   ::QuICC::SparseSM::Chebyshev::LinearMap::Y1 opY1(
      this->mspSetup->specSize() + 1, this->mspSetup->specSize(),
      this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.solver().setSpectralOperator(opY1.mat(), 1);
}

void D1Y1<base_t>::initBackend() const
{
   // Call parent initializer
   ILinearMapProjector::initBackend();

   // Initialize the solver
   this->mBackend.addSolver(1);
}

void D1Y1<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const Matrix>& in) const
{
   this->mBackend.input(tmp, in);
   auto specOp = this->mBackend.solver().getSpectralOperator();
   tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
   this->mBackend.getSolution(tmp, 1, 1);
}

void D1Y1<base_t>::applyPostOperator(Eigen::Ref<Matrix>) const {}

void D1Y1<base_t>::applyPreOperator(Matrix& tmp, const Eigen::Ref<const MatrixZ>& in,
   const bool useReal) const
{
   this->mBackend.input(tmp, in, useReal);
   auto specOp = this->mBackend.solver().getSpectralOperator();
   tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
   this->mBackend.getSolution(tmp, 1, 1);
}

void D1Y1<base_t>::applyPostOperator(Eigen::Ref<MatrixZ> rOut, const Matrix& tmp,
   const bool useReal) const
{
   this->mBackend.output(rOut, tmp, useReal);
}

} // namespace Projector
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
