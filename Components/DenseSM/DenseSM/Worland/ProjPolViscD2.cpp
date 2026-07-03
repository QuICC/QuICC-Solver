/**
 * @file ProjPolViscD2.cpp
 * @brief Implementation of the poloidal projection of the term -D2 u /rho
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/ProjPolViscD2.hpp"
#include "DenseSM/Worland/FC.hpp"
#include "DenseSM/Worland/DivR1D1FD1R1.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"

#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1dWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1WnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"


#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjPolViscD2::ProjPolViscD2(const int nNr, const int nNc, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{
   // bilaplacian
   this->mpOpA = std::make_shared<FC>(nNr, nNc, lOut, mOut, lF, mF,
         lIn, mIn, pF, alpha, dBeta);
   // F derivative term     
   this->mpOpB = std::make_shared<DivR1D1FD1R1>(nNr, nNc, lOut, mOut, lF, mF,
         lIn, mIn, pF, alpha, dBeta);
}

void ProjPolViscD2::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mpF->ls().size() != 1)
   {
      throw std::logic_error(
         "Operators are not implemented for forcing with multiple l");
   }

   const int lF = this->mpF->ls().at(0); // the real lf is actually 1 (r D2 / rho)

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
   Polynomial::Worland::slaplWnl slaplW;

   
   auto f = this->mpF->evaluate(igrid, lF, this->mMf); // I know this is actually r p(r^2); l=1, although lf=0

   Internal::Matrix opFFwdm1(igrid.size(), this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwdm1, opFFwdm1.cols(), lF+1, igrid, iweights,
      ev::Set());

   Internal::Matrix opFBwd(igrid.size(), this->mpF->nN());
   r_1W.compute<Internal::MHDFloat>(opFBwd, opFBwd.cols(),
      lF+1, igrid, Internal::Array(), ev::Set());

   Internal::Array r_1f = opFBwd * (opFFwdm1.transpose() * f); // l = 0 = lF

   // FIRST OPERATOR: (1/r)f* lapl(*) : l=lin+lf

   Internal::Matrix opLaplBwd(igrid.size(), this->cols());
   slaplW.compute<Internal::MHDFloat>(opLaplBwd, opLaplBwd.cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opPhys1 =  r_1f.asDiagonal() * opLaplBwd; // l = lin + lf = lin
   //Internal::Matrix opPhys1 =  opLaplBwd; // l = lin + lf = lin


   // SECOND OPERATOR: (1/r) D(f/r) D( r *) : l=lin+lf

   Internal::Matrix opFFwd(igrid.size(), this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, opFFwd.cols(), lF, igrid, iweights,
      ev::Set());

   Internal::Matrix opFDBwd(igrid.size(), this->mpF->nN());
   r_1dW.compute<Internal::MHDFloat>(opFDBwd, opFDBwd.cols(),
      lF, igrid, Internal::Array(), ev::Set());

   Internal::Array r_1Dr_1f =  opFDBwd * (opFFwd.transpose() * r_1f);  // l = 0 = lF

   Internal::Array spec1 = (opFFwd.transpose() * r_1f);

   Internal::Array spec2 = (opFFwd.transpose() * r_1Dr_1f);

   Internal::Matrix opBwdr_1drW(igrid.size(), this->cols());
   r_1drW.compute<Internal::MHDFloat>(opBwdr_1drW, opBwdr_1drW.cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());

   Internal::Matrix opPhys2 =  r_1Dr_1f.asDiagonal() * (igrid.asDiagonal() * opBwdr_1drW); // l = lin + lf = lin
   //Internal::Matrix opPhys2 =   (igrid.asDiagonal() * opBwdr_1drW); // l = lin + lf = lin

   // Combine the operators

   Internal::Matrix opFwdOut(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwdOut, opFwdOut.cols(), this->mLout, igrid,
      iweights, ev::Set());

   mat = opFwdOut.transpose() * (opPhys1 + opPhys2);
   //mat = -this->mpOpA->mat() + this->mpOpB->mat();
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
