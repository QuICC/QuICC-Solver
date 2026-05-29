/**
 * @file DivR2FD3R1.cpp
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
#include "DenseSM/Worland/DivR2FD3R1.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1WnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1dWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnlRecurrence.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

DivR2FD3R1::DivR2FD3R1(const int nNr, const int nNc, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{}

void DivR2FD3R1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3 * (2 * this->rows() + std::max(this->mLin, this->mLout) +
                          2 * this->mpF->nN() + this->mLf + 4)) /
                  2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::Wnl W;
   Polynomial::Worland::r_1Wnl<Polynomial::Worland::recurrence_t> r_1W;
   Polynomial::Worland::r_1drWnl<Polynomial::Worland::recurrence_t> r_1drW;
   Polynomial::Worland::r_1dWnl<Polynomial::Worland::recurrence_t> r_1dW;

   int lm;
   if(this->mLin==0)
   {
      lm = 1;
   }
   else
   {
      lm = this->mLin-1;
   }


   Internal::Matrix opBwd(igrid.size(), this->cols());
   r_1drW.compute<Internal::MHDFloat>(opBwd, this->cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());
   
   Internal::Matrix opFwd(igrid.size(), this->cols());
   W.compute<Internal::MHDFloat>(opFwd, this->cols(), this->mLin, igrid,
      iweights, ev::Set());
   Internal::Matrix opFwdlm(igrid.size(), this->cols());
   W.compute<Internal::MHDFloat>(opFwdlm, this->cols(), lm, igrid,
      iweights, ev::Set());

   Internal::Matrix opBwd1(igrid.size(), this->cols());
   r_1dW.compute<Internal::MHDFloat>(opBwd1, this->cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opBwd1lm(igrid.size(), this->cols());
   r_1dW.compute<Internal::MHDFloat>(opBwd1lm, this->cols(), lm, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opFFwd(igrid.size(), this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, this->mpF->nN(), this->mLf, igrid,
      iweights, ev::Set());
   Internal::Matrix opFBwd(igrid.size(), this->mpF->nN());
   r_1W.compute<Internal::MHDFloat>(opFBwd, this->mpF->nN(), this->mLf, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opFwdOut(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwdOut, this->rows(), this->mLout, igrid,
      iweights, ev::Set());

   //Internal::Array f = this->mpF->evaluate(igrid, this->mLf, this->mMf);
   Internal::Array f = opFBwd * opFFwd.transpose() * 
                       this->mpF->evaluate(igrid, this->mLf, this->mMf); // f/r in physical space

   mat = opFwdOut.transpose() * f.asDiagonal() 
            * opBwd1lm // r_1 d3(r*) in physical space
            * (opFwdlm.transpose() * igrid.asDiagonal()* opBwd1) // d2(r*) in spectral space
            * (opFwd.transpose() * igrid.asDiagonal() * opBwd); // d1(r*) in spectral space
   //copied from DivR2FD1R1
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
