/**
 * @file DivR1D1FD1R1C.cpp
 * @brief Source of the implementation of the spectral operator 1/r D(f
 * (-lapl(*)))
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/DivR1D1FD1R1C.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1dWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

DivR1D1FD1R1C::DivR1D1FD1R1C(const int nNr, const int nNc, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{}

void DivR1D1FD1R1C::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mpF->ls().size() != 1)
   {
      throw std::logic_error(
         "Operators are not implemented for forcing with multiple l");
   }

   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3 * (2 * this->rows() + std::max(this->mLin, this->mLout) +
                          2 * this->mpF->nN() + this->mLf + 4)) /
                  2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::slaplWnl slaplW;
   Polynomial::Worland::Wnl W;
   Polynomial::Worland::r_1dWnl<Polynomial::Worland::recurrence_t> r_1dW;
   Polynomial::Worland::r_1drWnl<Polynomial::Worland::recurrence_t> r_1drW;

   Internal::Matrix opBwd(igrid.size(), this->cols());
   slaplW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());
   opBwd = -opBwd;

   Internal::Matrix opFwdI(igrid.size(), this->cols());
   W.compute<Internal::MHDFloat>(opFwdI, opFwdI.cols(), this->mLin, igrid,
      iweights, ev::Set());

   Internal::Matrix opBwdD(igrid.size(), this->cols());
   r_1drW.compute<Internal::MHDFloat>(opBwdD, opBwdD.cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);

   Internal::Matrix opFFwd(igrid.size(), this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, opFFwd.cols(), this->mLf, igrid, iweights,
      ev::Set());

   Internal::Matrix opFBwd(igrid.size(), this->mpF->nN());
   r_1dW.compute<Internal::MHDFloat>(opFBwd, opFBwd.cols(),
      this->mLf, igrid, Internal::Array(), ev::Set());

   Internal::Array d1f = igrid.asDiagonal() * opFBwd * (opFFwd.transpose() * f);

   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, this->rows(), this->mLout, igrid,
      iweights, ev::Set());

   mat = opFwd.transpose() * d1f.asDiagonal() * opBwdD * opFwdI.transpose() *
          opBwd;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
