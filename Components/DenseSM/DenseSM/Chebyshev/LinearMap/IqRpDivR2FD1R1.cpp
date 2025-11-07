/**
 * @file IqRpDivR2FD1R1.cpp
 * @brief Source of the implementation of the spectral operator I^q r^p 1/I^q r^p f D(r
 * *)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR2FD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

IqRpDivR2FD1R1::IqRpDivR2FD1R1(const int nNr, const int nNc, const int q, const int p, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    IIqTripleHarmonicOperator(nNr, nNc, q, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP < 2)
   {
      throw std::logic_error("Radial prefactor needs to be at least I^q r^p");
   }

   if(this->mQ < 1)
   {
      throw std::logic_error("Quasi-inverse need to be at least 1");
   }
}

void IqRpDivR2FD1R1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   int rN = this->mpF->nN() + this->mP;

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();
   Matrix d1f = this->mpF->evaluateDiff(1, igrid, this->mLf, this->mMf, lb, ub);

   const Internal::Array& r = igrid;
   mat = Internal::Matrix::Zero(rows,cols);

   int fN = this->mpF->nN() + this->mP - 1;
   Matrix tf =
      (f.array()*r.array().pow(this->mP-1))
         .cast<MHDFloat>();
   Matrix cf = Utils::computeExpansion(tf, fN, lb, ub);

   Matrix fOp = Matrix::Zero(rows + fN + 2*this->mQ,cols);
   Utils::expansionProduct(fOp, fOp.rows(), fOp.cols(), cf, fN);
   mat = Utils::matIq(this->mQ, 1, rows, fOp.rows(), lb, ub) * fOp;

   if(this->mP > 2)
   {
      fN = this->mpF->nN() + this->mP - 2;
      tf =
         (-(r.array().pow(this->mP-2)*((this->mP-2)*f.array() + d1f.array()*r.array())))
            .cast<MHDFloat>();
   }
   else
   {
      fN = this->mpF->nN();
      tf =
         (-d1f.array()*r.array())
            .cast<MHDFloat>();
   }
   cf = Utils::computeExpansion(tf, fN, lb, ub);

   fOp.setZero();
   Utils::expansionProduct(fOp, fOp.rows(), fOp.cols(), cf, fN);
   mat += Utils::matIq(this->mQ, 0, rows, fOp.rows(), lb, ub) * fOp;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
