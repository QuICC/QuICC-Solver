/**
 * @file RpDivR1D1F.cpp
 * @brief Source of the implementation of the spectral operator r^p 1/r D(f *)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpDivR1D1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpDivR1D1F::RpDivR1D1F(const int nNr, const int nNc, const int p, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP != 4)
   {
      throw std::logic_error("Radial prefactor needs to be r^4");
   }
}

void RpDivR1D1F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   int fN = this->mpF->nN();
   int rN =
      2 * (std::max(this->rows(), this->cols()) + this->mpF->nN() + 4 + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();
   Matrix sF = Utils::computeExpansion(f, fN, lb, ub);

   Matrix tA = Matrix::Zero(rN,cols);
   Utils::expansionProduct(tA, cols+fN, cols, sF, fN);

   Matrix tB = Utils::evaluateD<1>(tA, this->cols() + fN, lb, ub);

   tB = igrid.array().pow(3).cast<MHDFloat>().matrix().asDiagonal() * tB;

   tA = Utils::computeExpansion(tB, this->rows(), lb, ub);
   mat = tA.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
