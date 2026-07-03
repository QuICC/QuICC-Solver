/**
 * @file ProjPolViscD1.cpp
 * @brief Implementation of the poloidal projection of the term D1 D1(u)/rho
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/ProjPolViscD1.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1dWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1WnlRecurrence.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjPolViscD1::ProjPolViscD1(const int nNr, const int nNc, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{
   /*
   // bilaplacian
   this->mpOpA = std::make_shared<FC2>(nNr, nNc, lOut, mOut, lF, mF,
         lIn, mIn, pF, alpha, dBeta);
   // F derivative term     
   this->mpOpB = std::make_shared<DivR1D1FD1R1C>(nNr, nNc, lOut, mOut, lF, mF,
         lIn, mIn, pF, alpha, dBeta);
   */
}

void ProjPolViscD1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mpF->ls().size() != 1)
   {
      throw std::logic_error(
         "Operators are not implemented for forcing with multiple l");
   }

   const int lF = this->mpF->ls().at(0);

   int lINm1;
   if (this->mLin>0)
   {
      lINm1 = this->mLin-1;
   }
   else if (this->mLin==0)
   {
      lINm1 = 1;
   }
   else
   {
      throw std::logic_error(
         "Attempted Worland transform with lIn<0");
   }

   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3 * (2 * this->rows() + std::max(this->mLin, this->mLout) +
                          2 * this->mpF->nN() + this->mLf + 4)) /2;

   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::Wnl W;
   Polynomial::Worland::r_1dWnl<Polynomial::Worland::recurrence_t> r_1dW;
   Polynomial::Worland::r_1drWnl<Polynomial::Worland::recurrence_t> r_1drW;
   Polynomial::Worland::r_1Wnl<Polynomial::Worland::recurrence_t> r_1W;

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf); // I know this is actually r^2 p(r^2); l=2, although lf=0

   // FIRST OPERATOR: f* D((1/r)*) : l=lin+lf-2

   Internal::Matrix opBwd_r_1W(igrid.size(), this->cols());
   r_1W.compute<Internal::MHDFloat>(opBwd_r_1W, opBwd_r_1W.cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opFwdm1(igrid.size(), this->cols());
   W.compute<Internal::MHDFloat>(opFwdm1, opFwdm1.cols(), lINm1, igrid,
      iweights, ev::Set());

   Internal::Matrix opBwd_r_1dW(igrid.size(), this->cols());
   r_1dW.compute<Internal::MHDFloat>(opBwd_r_1dW, opBwd_r_1dW.cols(), lINm1, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opFFwd(igrid.size(), this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, opFFwd.cols(), lF, igrid, iweights,
      ev::Set());

   Internal::Matrix opFBwd_r_1(igrid.size(), this->mpF->nN());
   r_1W.compute<Internal::MHDFloat>(opFBwd_r_1, opFBwd_r_1.cols(), lF, igrid,
      Internal::Array(), ev::Set());

   Internal::Array rf = igrid.asDiagonal() * f;

   // returns an l = lin -2 + lf. I know f is actually l=2. RETURN l=lin=lin+lf
   Internal::Matrix opPhys1 = rf.asDiagonal() * // rf
                              opBwd_r_1dW* opFwdm1.transpose() * opBwd_r_1W; // (1/r)D(1/r *)


   // SECOND OPERATOR: r D( f D( (1/r)D(r*) ) ): l=lin+lf-2
   Internal::Matrix opBwdr_r_1drW(igrid.size(), this->cols());
   r_1drW.compute<Internal::MHDFloat>(opBwdr_r_1drW, opBwdr_r_1drW.cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opBwdr_r_1dW(igrid.size(), this->cols());
   r_1dW.compute<Internal::MHDFloat>(opBwdr_r_1dW, opBwdr_r_1dW.cols(), lINm1, igrid,
      Internal::Array(), ev::Set());

   //supposedly this takes an l = lin -3 + lF. 
   // But I know I am passing something that is actually lf=2. so I can write lin-1+lf (lf in input being zero)
   Internal::Matrix opFFwd1(igrid.size(), this->cols()+this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd1, opFFwd1.cols(), lINm1 + lF, igrid, iweights,
      ev::Set());

   Internal::Matrix opFBwdr_r_1drW(igrid.size(), this->cols()+this->mpF->nN());
   r_1drW.compute<Internal::MHDFloat>(opFBwdr_r_1drW, opFBwdr_r_1drW.cols(), lINm1 + lF, igrid,
      Internal::Array(), ev::Set());
   
   Internal::Matrix opPhys2 = igrid.asDiagonal() * (igrid.asDiagonal() * // r^2
                              opFBwdr_r_1drW * opFFwd1.transpose() * // (1/r)D(r *)
                              f.asDiagonal() * // f      //(f.array()).matrix().asDiagonal() * // r*f igrid.array() * 
                              opBwdr_r_1dW * opFwdm1.transpose() * opBwdr_r_1drW); // (1/r)d( (1/r)d(r *) )

   // OUTPUT OPERATOR
   // I know  opPhys1 and  opPhys2 are lIn individually, but 
   // -Lout2 * opPhys1 + opPhys2 has to be lIn+2 (or say, lIn+2 + lF)

   Internal::Matrix opFFwd2(igrid.size(), this->cols()+this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd2, opFFwd2.cols(), this->mLin + lF + 2, igrid, iweights,
      ev::Set());

   Internal::Matrix opBwdOut2(igrid.size(), this->cols()+this->mpF->nN());
   r_1W.compute<Internal::MHDFloat>(opBwdOut2, opBwdOut2.cols(), this->mLin + lF + 2, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opFFwd3(igrid.size(), this->cols()+this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd3, opFFwd3.cols(), this->mLin + lF + 1, igrid, iweights,
      ev::Set());

   Internal::Matrix opBwdOut1(igrid.size(), this->cols()+this->mpF->nN());
   r_1W.compute<Internal::MHDFloat>(opBwdOut1, opBwdOut1.cols(), this->mLin + lF + 1, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opFwdOut(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwdOut, opFwdOut.cols(), this->mLout, igrid,
      iweights, ev::Set());
   
   const int Lout2 = this->mLout*(this->mLout+1);

   mat = opFwdOut.transpose() * opBwdOut1 * opFFwd3.transpose() * opBwdOut2 * opFFwd2.transpose() * (-Lout2 * opPhys1 + opPhys2);
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
