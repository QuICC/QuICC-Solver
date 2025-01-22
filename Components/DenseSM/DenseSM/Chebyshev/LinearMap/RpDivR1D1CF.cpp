/**
 * @file RpDivR1D1CF.cpp
 * @brief Source of the implementation of the spectral operator r^p 1/r
 * D(-lapl(f) *) multiplied by gaunt coefficient
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpDivR1D1CF.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpDivR1D1CF::RpDivR1D1CF(const int nNr, const int nNc, const int p, const int lOut,
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

void RpDivR1D1CF::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   const auto pId = GridPurpose::SIMULATION;
   int rN =
      2 * (std::max(this->rows(), this->cols()) + this->mpF->nN() + 4 + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix tA = Matrix::Identity(rN, this->cols());

   Matrix tB = Utils::evaluate(tA, this->cols(), lb, ub);
   Matrix td1B = Utils::evaluateD<1>(tA, this->cols(), lb, ub);
   Matrix td2B = Utils::evaluateD<2>(tA, this->cols(), lb, ub);
   Matrix td3B = Utils::evaluateD<3>(tA, this->cols(), lb, ub);

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();
   Matrix d1f = this->mpF->evaluateDiff(1, igrid, this->mLf, this->mMf, lb, ub);
   Matrix d2f = this->mpF->evaluateDiff(2, igrid, this->mLf, this->mMf, lb, ub);
   Matrix d3f = this->mpF->evaluateDiff(3, igrid, this->mLf, this->mMf, lb, ub);

   const int l = this->mLf;
   const Internal::Array& ir = igrid;
   Array r = ir.cast<MHDFloat>();

   tA =
      (-f).asDiagonal() * l * (1 + l) * ((-r).asDiagonal() * td1B + 2.0 * tB) +
      r.asDiagonal() *
         (((-r).array() * (2.0 * d1f.array() + r.array() * d2f.array()))
                  .matrix()
                  .asDiagonal() *
               td1B +
            (d1f.array() * (2 + l + l * l) -
               r.array() * (2.0 * d2f.array() + d3f.array() * r.array()))
                  .matrix()
                  .asDiagonal() *
               tB);

   Matrix tC = Utils::computeExpansion(tA, this->rows(), lb, ub);

   mat = tC.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
