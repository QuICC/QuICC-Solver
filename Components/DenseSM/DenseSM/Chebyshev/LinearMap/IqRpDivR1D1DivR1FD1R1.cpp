/**
 * @file IqRpDivR1D1DivR1FD1R1.cpp
 * @brief Source of the implementation of the spectral operator I^q r^p 1/r D(1/r
 * D(r f *))
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR1D1DivR1FD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

IqRpDivR1D1DivR1FD1R1::IqRpDivR1D1DivR1FD1R1(const int nNr, const int nNc, const int q, const int p,
   const int lOut, const int mOut, const int lF, const int mF, const int lIn,
   const int mIn, std::shared_ptr<RadialTorPolFunction> pF,
   const Scalar_t lower, const Scalar_t upper) :
    IIqTripleHarmonicOperator(nNr, nNc, q, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP != 4)
   {
      throw std::logic_error("Radial prefactor needs to be r^4");
   }
}

void IqRpDivR1D1DivR1FD1R1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   int rN = this->mpF->nN() + 4;

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();
   Matrix d1f = this->mpF->evaluateDiff(1, igrid, this->mLf, this->mMf, lb, ub);

   const Internal::Array& r = igrid;
   mat = Internal::Matrix::Zero(rows,cols);

   int fN = this->mpF->nN() + 3;
   Matrix tf =
      (f.array()*r.array().pow(3))
         .cast<MHDFloat>();
   Matrix cf = Utils::computeExpansion(tf, fN, lb, ub);

   Matrix fOp = Matrix::Zero(rows + fN + 2*this->mQ,cols);
   Utils::expansionProduct(fOp, fOp.rows(), fOp.cols(), cf, fN);
   mat = Utils::matIq(this->mQ, 2, rows, fOp.rows(), lb, ub) * fOp;

   fN = this->mpF->nN() + 2;
   tf =
      (-(r.array().pow(2)*(5.0*f.array() + d1f.array()*r.array())))
         .cast<MHDFloat>();
   cf = Utils::computeExpansion(tf, fN, lb, ub);

   fOp.setZero();
   Utils::expansionProduct(fOp, fOp.rows(), fOp.cols(), cf, fN);
   mat += Utils::matIq(this->mQ, 1, rows, fOp.rows(), lb, ub) * fOp;

   fN = this->mpF->nN() + 1;
   tf =
      (3.0*r.array()*(f.array() + d1f.array()* r.array()))
         .cast<MHDFloat>();
   cf = Utils::computeExpansion(tf, fN, lb, ub);

   fOp.setZero();
   Utils::expansionProduct(fOp, fOp.rows(), fOp.cols(), cf, fN);
   mat += Utils::matIq(this->mQ, 0, rows, fOp.rows(), lb, ub) * fOp;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
