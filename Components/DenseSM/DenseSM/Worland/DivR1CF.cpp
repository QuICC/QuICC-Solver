/**
 * @file DivR1CF.cpp
 * @brief Source of the implementation of the spectral operator -lapl(f)/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/DivR1CF.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1WnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

DivR1CF::DivR1CF(const int nNr, const int nNc, const int lOut, const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha, const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha, dBeta)
{}

void DivR1CF::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if(this->mpF->ls().size() != 1)
   {
      throw std::logic_error("Operators are not implemented for forcing with multiple l");
   }
   const int lF = this->mpF->ls().at(0);

   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3*(2*(this->rows() + this->mpF->nN()) + std::max(this->mLin, this->mLout)))/2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);
   
   Polynomial::Worland::r_1Wnl<Polynomial::Worland::recurrence_t> r_1W;
   Internal::Matrix opBwd(igrid.size(), this->cols());
   r_1W.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mLin, igrid, Internal::Array(), ev::Set());
   
   Polynomial::Worland::Wnl W;
   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mLout, igrid, iweights, ev::Set());

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);
   Internal::Matrix opFFwd(igrid.size(), this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, opFFwd.cols(), lF, igrid, iweights, ev::Set());
   Polynomial::Worland::slaplWnl slaplW;
   Internal::Matrix opFBwd(igrid.size(), this->mpF->nN());
   slaplW.compute<Internal::MHDFloat>(opFBwd, opFBwd.cols(), lF, igrid, Internal::Array(), ev::Set());
   Array cf = -(opFBwd * (opFFwd.transpose() * f));

   mat = opFwd.transpose() * cf.asDiagonal() * opBwd;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
