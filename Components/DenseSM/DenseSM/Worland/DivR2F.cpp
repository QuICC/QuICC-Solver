/**
 * @file DivR2F.cpp
 * @brief Source of the implementation of the spectral operator f/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/DivR2F.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1WnlRecurrence.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

DivR2F::DivR2F(const int nNr, const int nNc, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{}

void DivR2F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3 * (2 * this->rows() + std::max(this->mLin, this->mLout) +
                          2 * this->mpF->nN() + this->mLf + 4)) /
                  2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::r_1Wnl<Polynomial::Worland::recurrence_t> r_1wnl;
   Internal::Matrix opBwd(igrid.size(), this->cols());
   r_1wnl.compute<Internal::MHDFloat>(opBwd, this->cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());

   Polynomial::Worland::Wnl wnl;
   Internal::Matrix opFFwd(igrid.size(), this->mpF->nN());
   wnl.compute<Internal::MHDFloat>(opFFwd, this->mpF->nN(), this->mLf, igrid,
      iweights, ev::Set());
   Internal::Matrix opFBwd(igrid.size(), this->mpF->nN());
   r_1wnl.compute<Internal::MHDFloat>(opFBwd, this->mpF->nN(), this->mLf, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opFwd(igrid.size(), this->rows());
   wnl.compute<Internal::MHDFloat>(opFwd, this->rows(), this->mLout, igrid,
      iweights, ev::Set());

   Internal::Array f = opFBwd * opFFwd.transpose() *
                       this->mpF->evaluate(igrid, this->mLf, this->mMf);

   mat = opFwd.transpose() * f.asDiagonal() * opBwd;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
