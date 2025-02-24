/**
 * @file R0F.cpp
 * @brief Source of the implementation of the spectral operator f * 
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R0F.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R0F::R0F(const int nNr, const int nNc, const int p, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   /*
   if(this->mP < 2)
   {
      throw std::logic_error("Radial prefactor should be at least r^p");
   }*/
}

void R0F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   int fN = this->mpF->nN()- 1; //+ this->mP;
   int rN =
      2 * (fN + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);
   //f = igrid.array().pow(this->mP-1).matrix().asDiagonal() * f;
   Matrix gF = f.cast<MHDFloat>();

   Matrix sF = Utils::computeExpansion(gF, fN, lb, ub);

   mat = Matrix::Zero(rows,cols);
   Utils::expansionProduct(mat, rows, cols, sF, fN);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
