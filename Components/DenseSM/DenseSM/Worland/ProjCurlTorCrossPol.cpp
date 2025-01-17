/**
 * @file ProjCurlTorCrossPol.cpp
 * @brief Source of the implementation of the projection r Curl Tor ^ Tor
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/ProjCurlTorCrossPol.hpp"
#include "DenseSM/Worland/DivR1F.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlTorCrossPol::ProjCurlTorCrossPol(const int nNr, const int nNc, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB,
   std::shared_ptr<RadialTorPolFunction> pTorA, std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t alpha, const Scalar_t dBeta) :
    IProjCrossOperator(nNr, nNc, lOut, mOut, lA, mA, lB, mB, alpha, dBeta)
{
   // Radial function A is given
   if(pTorA && pPolB == nullptr)
   {
      this->mpOp = std::make_shared<DivR1F>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pTorA, alpha, dBeta);
   }
   // Radial function B is given
   else if(pPolB && pTorA == nullptr)
   {
      this->mpOp = std::make_shared<DivR1F>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, alpha, dBeta);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }

   this->mIsZero = (this->elsasser(this->mLa, this->mMa, this->mLb, this->mMb, this->mLout, this->mMout) == 0);
}

void ProjCurlTorCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if(this->mIsZero)
   {
      mat = Matrix::Zero(this->rows(), this->cols());
   }
   else
   {
      auto&& lg = this->mLout;
      auto&& la = this->mLa;
      auto&& lb = this->mLb;

      const MHDFloat L2b = static_cast<MHDFloat>(lb*(lb + 1));

      const MHDFloat Labg = this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

      MHDFloat c = std::pow(-1.0, la  + lb + lg - 1)*L2b*Labg;

      mat = c*this->mpOp->mat();
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
