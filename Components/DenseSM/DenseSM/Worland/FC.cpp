/**
 * @file FC.cpp
 * @brief Source of the implementation of the spectral operator f(-lapl(*))
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/FC.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Worland/Tags.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

FC::FC(const int nNr, const int nNc, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, alpha,
       dBeta)
{}

void FC::buildOpImpl(Internal::Matrix& mat, const int rows,
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
                          2 * this->mpF->nN() + this->mLf)) /
                  2;
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, nR);

   Polynomial::Worland::slaplWnl slaplW;
   Internal::Matrix opBwd(igrid.size(), this->cols());
   slaplW.compute<Internal::MHDFloat>(opBwd, opBwd.cols(), this->mLin, igrid,
      Internal::Array(), ev::Set());
   opBwd = -opBwd;

   Polynomial::Worland::Wnl W;
   Internal::Matrix opFwd(igrid.size(), this->rows());
   W.compute<Internal::MHDFloat>(opFwd, opFwd.cols(), this->mLout, igrid,
      iweights, ev::Set());

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);

   mat = opFwd.transpose() * f.asDiagonal() * opBwd;
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
