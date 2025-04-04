/**
 * @file EnergySLaplR2.cpp
 * @brief Source of the implementation of the Chebyshev energy r^2 spherical
 * laplacian reductor, with linear map y = ax + b
 */

// System includes
//
#include <cassert>

// Project includes
//
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Reductor/Base/EnergySLaplR2.hpp"
#include "QuICC/Polynomial/Quadrature/ChebyshevRule.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/I1.hpp"
#include "QuICC/SparseSM/Chebyshev/LinearMap/Y2.hpp"

#include <iostream>
namespace QuICC {

namespace Transform {

namespace Fft {

namespace Chebyshev {

namespace LinearMap {

namespace Reductor {

void EnergySLaplR2<base_t>::initOperator() const
{
   // Check for division by 0!
   assert(this->mspSetup->lower() > 0.0 || this->mspSetup->upper() < 0.0);

   ::QuICC::SparseSM::Chebyshev::LinearMap::I1 spasmI1(
      this->mspSetup->specSize() + 3, this->mspSetup->specSize() + 3,
      this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.solver().setOperator(spasmI1.mat(), 2, 2);

   ::QuICC::SparseSM::Chebyshev::LinearMap::Y2 spasmY2(
      this->mspSetup->specSize() + 2, this->mspSetup->specSize() + 2,
      this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.solver().setSpectralOperator(spasmY2.mat(), 1);

   Internal::Array igrid, iweights;
   Polynomial::Quadrature::ChebyshevRule quad;
   quad.computeQuadrature(igrid, iweights, 2*this->mspSetup->fwdSize(),
      this->mspSetup->lower(), this->mspSetup->upper());
   this->mBackend.setScaler(igrid.array().pow(-1).cast<MHDFloat>().matrix());
}

void EnergySLaplR2<base_t>::initBackend() const
{
   // Call parent initializer
   ILinearMapEnergy::initBackend();

   // Initialize the solver
   this->mBackend.addSolver(2);
}

void EnergySLaplR2<base_t>::applyPreOperator(Matrix& tmp, const Matrix& in) const
{
   Matrix tmp2(tmp.rows(), tmp.cols());
   this->mBackend.input(tmp2, in);

   this->mBackend.input(tmp, in, 1);
   this->mBackend.getSolution(tmp, 3, -1);
   auto specOp = this->mBackend.solver().getSpectralOperator();
   tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
   this->mBackend.getSolution(tmp, 1, 2);

   int col = 0;
   int mult = 0;
   for(int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      MHDFloat l = static_cast<MHDFloat>(this->mspSetup->slow(i));
      MHDFloat ll1 = l * (l + 1.0);
      mult = this->mspSetup->mult(i);
      tmp.block(0, col, this->mspSetup->specSize(), mult) -= ll1*tmp2.block(0, col, this->mspSetup->specSize(), mult);
      col += mult;
   }
   assert(col == this->mspSetup->blockSize());
}

void EnergySLaplR2<base_t>::applyPostOperator(Matrix& rOut, const Matrix& tmp) const
{
   assert(rOut.cols() == 1);
   this->mBackend.output(rOut, tmp);
}

void EnergySLaplR2<base_t>::applyPreOperator(Matrix& tmp, const MatrixZ& in,
   const bool useReal) const
{
   Matrix tmp2(tmp.rows(), tmp.cols());
   this->mBackend.input(tmp2, in, useReal);

   this->mBackend.input(tmp, in, 1, useReal);
   this->mBackend.getSolution(tmp, 3, -1);
   auto specOp = this->mBackend.solver().getSpectralOperator();
   tmp.topRows(specOp.rows()) = specOp * tmp.topRows(specOp.cols());
   this->mBackend.getSolution(tmp, 1, 2);

   int col = 0;
   int mult = 0;
   for(int i = 0; i < this->mspSetup->slowSize(); i++)
   {
      MHDFloat l = static_cast<MHDFloat>(this->mspSetup->slow(i));
      MHDFloat ll1 = l * (l + 1.0);
      mult = this->mspSetup->mult(i);
      tmp.block(0, col, this->mspSetup->specSize(), mult) -= ll1*tmp2.block(0, col, this->mspSetup->specSize(), mult);
      col += mult;
   }
   assert(col == this->mspSetup->blockSize());
}

} // namespace Reductor
} // namespace LinearMap
} // namespace Chebyshev
} // namespace Fft
} // namespace Transform
} // namespace QuICC
