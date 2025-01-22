/**
 * @file RpDivR1CF.cpp
 * @brief Source of the implementation of the spectral operator r^p 1/r
 * (-lapl(f)) *
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpDivR1CF.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpDivR1CF::RpDivR1CF(const int nNr, const int nNc, const int p, const int lOut,
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

void RpDivR1CF::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   const auto pId = GridPurpose::SIMULATION;
   int fN = this->mpF->nN() + 1;
   int rN =
      2 * (std::max(this->rows(), this->cols()) + this->mpF->nN() + 4 + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();
   Matrix d1f = this->mpF->evaluateDiff(1, igrid, this->mLf, this->mMf, lb, ub);
   Matrix d2f = this->mpF->evaluateDiff(2, igrid, this->mLf, this->mMf, lb, ub);
   Matrix d3f = this->mpF->evaluateDiff(3, igrid, this->mLf, this->mMf, lb, ub);

   const int l = this->mLf;
   const Internal::Array& r = igrid;

   f =
      (r.array() * (f.array() * l * (1 + l) -
                      r.array() * (2 * d1f.array() + d2f.array() * r.array())))
         .cast<MHDFloat>();
   d1f = Utils::computeExpansion(f, fN, lb, ub);

   mat = Internal::Matrix::Zero(rows,cols);
   Utils::expansionProduct(mat, mat.rows(), mat.cols(), d1f, fN);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
