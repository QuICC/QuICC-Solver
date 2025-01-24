/**
 * @file IqRpDivR1F.cpp
 * @brief Source of the implementation of the spectral operator I^q r^p f/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/Typedefs.hpp"

#include <iostream>
namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

IqRpDivR1F::IqRpDivR1F(const int nNr, const int nNc, const int q, const int p, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    IIqTripleHarmonicOperator(nNr, nNc, q, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP < 2)
   {
      throw std::logic_error("Radial prefactor should be at least I^q r^p");
   }

   if(this->mQ < 1)
   {
      throw std::logic_error("Quasi-inverse need to be at least 1");
   }
}

void IqRpDivR1F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   int rN = this->mpF->nN() + this->mP;

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Internal::Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf);

   const Internal::Array& r = igrid;
   mat = Internal::Matrix::Zero(rows,cols);

   int fN = this->mpF->nN() + this->mP-1;
   Matrix tf =
      (f.array()*r.array().pow(this->mP-1))
         .cast<MHDFloat>();
   Matrix cf = Utils::computeExpansion(tf, fN, lb, ub);
   std::cerr << "fN: " << fN << std::endl;
   std::cerr << cf << std::endl;

   Matrix fOp = Matrix::Zero(rows + fN + 2*this->mQ,cols);
   Utils::expansionProduct(fOp, fOp.rows(), fOp.cols(), cf, fN);
   std::cerr << fOp << std::endl;
   std::cerr << Utils::matIq(this->mQ, 0, rows, fOp.rows(), lb, ub) << std::endl;
   mat = Utils::matIq(this->mQ, 0, rows, fOp.rows(), lb, ub) * fOp;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
