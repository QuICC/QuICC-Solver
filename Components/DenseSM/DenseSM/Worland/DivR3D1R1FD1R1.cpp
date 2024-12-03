/**
 * @file DivR3D1R1FD1R1.hpp.cpp
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
#include "DenseSM/Worland/DivR3D1R1FD1R1.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/r_1WnlRecurrence.hpp"
#include "QuICC/Polynomial/Worland/r_1drWnlRecurrence.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

#include<iostream>
namespace QuICC {

namespace DenseSM {

namespace Worland {

DivR3D1R1FD1R1::DivR3D1R1FD1R1(const int nNr, const int nNc, const int lOut, const int lF, const int lIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha, const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, lF, lIn, pF, alpha, dBeta)
{}

void DivR3D1R1FD1R1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   std::cerr << "lOut : " << this->mLout << std::endl;
   std::cerr << "lF : " << this->mLf << std::endl;
   std::cerr << "lIn : " << this->mLin << std::endl;
   std::cerr << "F nN : " << this->mpF->nN() << std::endl;
   std::cerr << "nNr : " << this->rows() << std::endl;
   std::cerr << "nNc : " << this->cols() << std::endl;
   namespace ev = Polynomial::Worland::Evaluator;
   const int nR = (3*(2*this->rows() + std::max(this->mLin, this->mLout)))/2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);
   
   Polynomial::Worland::Wnl W;
   Polynomial::Worland::r_1Wnl<Polynomial::Worland::recurrence_t> r_1W;
   Polynomial::Worland::r_1drWnl<Polynomial::Worland::recurrence_t> r_1drW;

   Internal::Matrix opBwd(igrid.size(), this->cols());
   r_1drW.compute<Internal::MHDFloat>(opBwd, this->cols(), this->mLin, igrid, Internal::Array(), ev::Set());

   Internal::Matrix opFFwd(igrid.size(), this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFwd, opFFwd.cols(), this->mLf, igrid, iweights, ev::Set());
   Internal::Matrix opFBwd(igrid.size(), this->mpF->nN());
   r_1drW.compute<Internal::MHDFloat>(opFBwd, opFBwd.cols(), this->mLf, igrid, Internal::Array(), ev::Set());

   Internal::Matrix opFFFwd(igrid.size(), this->cols() + this->mpF->nN());
   W.compute<Internal::MHDFloat>(opFFFwd, opFFFwd.cols(), this->mLin + this->mLf - 2, igrid, iweights, ev::Set());
   Internal::Matrix opFFBwd(igrid.size(), this->cols() + this->mpF->nN());
   r_1W.compute<Internal::MHDFloat>(opFFBwd, opFFBwd.cols(), this->mLin + this->mLf - 2, igrid, Internal::Array(), ev::Set());
   
   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, this->rows(), this->mLout, igrid, iweights, ev::Set());

   Internal::Array f = opFBwd * opFFwd.transpose() * this->mpF->evaluate(igrid, this->mLf);

   mat = opFwd.transpose() * opFFBwd * opFFFwd.transpose() * (f.asDiagonal() * opBwd);
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
