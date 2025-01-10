/**
 * @file DivR2D1R1FC.hpp.cpp
 * @brief Source of the implementation of the spectral operator 1/r D(r f) (-lapl(*))/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/DivR2D1R1FC.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1WnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

DivR2D1R1FC::DivR2D1R1FC(const int nNr, const int nNc, const int lOut, const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha, const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha, dBeta)
{}

void DivR2D1R1FC::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3*(2*this->rows() + std::max(this->mLin, this->mLout) + 2*this->mpF->nN() + this->mLf + 4))/2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::Wnl W;
   Polynomial::Worland::slaplWnl slaplW;
   Polynomial::Worland::r_1Wnl<Polynomial::Worland::recurrence_t> r_1W;
   Polynomial::Worland::r_1drWnl<Polynomial::Worland::recurrence_t> r_1drW;
   
   Internal::Matrix opCIBwd(igrid.size(), this->cols());
   slaplW.compute<Internal::MHDFloat>(opCIBwd, opCIBwd.cols(), this->mLin, igrid, Internal::Array(), ev::Set());
   opCIBwd = -opCIBwd;

   Internal::Matrix opIFwd(igrid.size(), this->cols());
   W.compute<Internal::MHDFloat>(opIFwd, opIFwd.cols(), this->mLin, igrid, iweights, ev::Set());
   
   Internal::Matrix opBwd(igrid.size(), this->cols());
   r_1W.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mLin, igrid, Internal::Array(), ev::Set());
   
   Internal::Matrix opFFwd(igrid.size(), this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, opFFwd.cols(), this->mLf, igrid, iweights, ev::Set());
   Internal::Matrix opFBwd(igrid.size(), this->mpF->nN());
   r_1drW.compute<Internal::MHDFloat>(opFBwd, opFBwd.cols(), this->mLf, igrid, Internal::Array(), ev::Set());

   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mLout, igrid, iweights, ev::Set());

   Internal::Array f = opFBwd * opFFwd.transpose() * this->mpF->evaluate(igrid, this->mLf, this->mMf);

   mat = opFwd.transpose() * f.asDiagonal() * opBwd * (opIFwd.transpose() * opCIBwd);
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
