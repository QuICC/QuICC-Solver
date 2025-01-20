/**
 * @file DivR1D1FC.hpp.cpp
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
#include "DenseSM/Worland/DivR1D1FC.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1dWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

DivR1D1FC::DivR1D1FC(const int nNr, const int nNc, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{}

void DivR1D1FC::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mpF->ls().size() != 1)
   {
      throw std::logic_error(
         "Operators are not implemented for forcing with multiple l");
   }
   const int lF = this->mpF->ls().at(0);

   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3 * (2 * this->rows() + std::max(this->mLin, this->mLout) +
                          2 * this->mpF->nN() + this->mLf + 4)) /
                  2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::slaplWnl slaplW;
   Internal::Matrix opBwd(igrid.size(), this->cols());
   slaplW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());
   opBwd = -opBwd;

   Polynomial::Worland::Wnl W;
   Polynomial::Worland::r_1dWnl<Polynomial::Worland::recurrence_t> r_1dW;

   Internal::Matrix opFFwd(igrid.size(), this->cols() + this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, opFFwd.cols(), this->mLin + this->mLf,
      igrid, iweights, ev::Set());
   Internal::Matrix opFBwd(igrid.size(), this->cols() + this->mpF->nN());
   r_1dW.compute<Internal::MHDFloat>(opFBwd, opFBwd.cols(),
      this->mLin + this->mLf, igrid, Internal::Array(), ev::Set());

   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, this->rows(), this->mLout, igrid,
      iweights, ev::Set());

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);

   mat = opFwd.transpose() * opFBwd * opFFwd.transpose() *
         (f.asDiagonal() * opBwd);
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
