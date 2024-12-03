/**
 * @file DivR1D1DivR1FD1R1.hpp.cpp
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
#include "DenseSM/Worland/DivR1D1DivR1FD1R1.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1dWnlRecurrence.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

DivR1D1DivR1FD1R1::DivR1D1DivR1FD1R1(const int nNr, const int nNc, const int lOut, const int lF, const int lIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha, const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, lF, lIn, pF, alpha, dBeta)
{}

void DivR1D1DivR1FD1R1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3*(2*this->rows() + std::max(this->mLin, this->mLout)))/2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);
   
   Polynomial::Worland::r_1drWnl<Polynomial::Worland::recurrence_t> r_1drW;
   Polynomial::Worland::r_1dWnl<Polynomial::Worland::recurrence_t> r_1dW;
   Polynomial::Worland::Wnl W;

   Internal::Matrix opBwd(igrid.size(), this->cols());
   r_1drW.compute<Internal::MHDFloat>(opBwd, this->cols(), this->mLin, igrid, Internal::Array(), ev::Set());

   Internal::Matrix opFFwd(igrid.size(), this->cols() + this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, opFFwd.cols(), this->mLf + this->mLf - 1, igrid, iweights, ev::Set());
   Internal::Matrix opFBwd(igrid.size(), this->cols() + this->mpF->nN());
   r_1dW.compute<Internal::MHDFloat>(opFBwd, opFBwd.cols(), this->mLf + this->mLf - 1, igrid, Internal::Array(), ev::Set());

   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, this->rows(), this->mLout, igrid, iweights, ev::Set());

   Internal::Array f = this->mpF->evaluate(igrid, this->mLf);

   mat = opFwd.transpose() * opFBwd * opFFwd.transpose() * f.asDiagonal() * opBwd;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
