/**
 * @file RpDivR2D1R1FC.cpp
 * @brief Source of the implementation of the spectral operator r^p 1/r D(r f)
 * (-lapl(*))/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpDivR2D1R1FC.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpDivR2D1R1FC::RpDivR2D1R1FC(const int nNr, const int nNc, const int p, const int lOut,
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

void RpDivR2D1R1FC::buildOpImpl(Internal::Matrix& mat, const int rows,
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

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();
   Matrix d1f = this->mpF->evaluateDiff(1, igrid, this->mLf, this->mMf, lb, ub);

   const int l = this->mLin;
   const Internal::Array& r = igrid;

   tA = -((f.array() + d1f.array() * r.cast<MHDFloat>().array())
             .matrix()
             .asDiagonal() *
          (r.cast<MHDFloat>().asDiagonal() *
                (2.0 * td1B + r.cast<MHDFloat>().asDiagonal() * td2B) -
             l * (1.0 + l) * tB));

   tB = Utils::computeExpansion(tA, this->rows(), lb, ub);
   mat = tB.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
