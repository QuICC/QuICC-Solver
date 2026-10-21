/**
 * @file SLapl.cpp
 * @brief Source of the implementation of spherical laplacian tau operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/SLapl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

SLapl::SLapl(const int rows, const int cols, const Scalar_t alpha,
   const Scalar_t dBeta, const int l) :
    ITauOperator(rows, cols, alpha, dBeta, l)
{
}

void SLapl::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (2 * this->rows() + this->mL);

   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::slaplWnl bW;
   Internal::Matrix opBwd(igrid.size(), this->rows());
   bW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mL, igrid,
         Internal::Array(), ev::Set());

   Polynomial::Worland::Wnl W;
   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mL, igrid,
         iweights, ev::Set());

   mat = opFwd.transpose() * opBwd;

   // Prune zeros
   for(int i = 0; i < mat.cols(); i++)
   {
      mat.block(i, i, mat.rows() - i, 1).setZero();
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
