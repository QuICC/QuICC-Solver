/**
 * @file RpDivR2D1R1F.cpp
 * @brief Source of the implementation of the spectral operator r^p 1/r^p D(r f)
 * *
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpDivR2D1R1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpDivR2D1R1F::RpDivR2D1R1F(const int nNr, const int nNc, const int p, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP < 2)
   {
      throw std::logic_error("Radial prefactor needs to be at least r^p");
   }
}

void RpDivR2D1R1F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;

   int fN = this->mpF->nN() + this->mP - 2;
   int rN =
      2 * (fN + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix fA =
      this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();
   Matrix fB = Utils::computeExpansion(fA, fN, lb, ub);

   fA = Utils::evaluateOp<cheb::Projector::D1Y1<cheb::base_t>>(fB, this->mpF->nN(), lb, ub);
   if(this->mP - 2 > 0)
   {
      fA = igrid.array().pow(this->mP-2).cast<MHDFloat>().matrix().asDiagonal() * fA;
   }
   fB = Utils::computeExpansion(fA, fN, lb, ub);

   mat = Matrix::Zero(rows,cols);
   Utils::expansionProduct(mat, rows, cols, fB, fN);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
