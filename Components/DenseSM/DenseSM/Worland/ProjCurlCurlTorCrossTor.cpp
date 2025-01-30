/**
 * @file ProjCurlCurlTorCrossTor.cpp
 * @brief Source of the implementation of the projection r Curl Curl Tor ^ Tor
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/DivR1F.hpp"
#include "DenseSM/Worland/ProjCurlCurlTorCrossTor.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlCurlTorCrossTor::ProjCurlCurlTorCrossTor(const int nNr, const int nNc,
   const int q, const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, std::shared_ptr<RadialTorPolFunction> pTorA,
   std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t alpha,
   const Scalar_t dBeta) :
    IProjCrossOperator(nNr, nNc, q, lOut, mOut, lA, mA, lB, mB, alpha, dBeta)
{
   auto prefactor = [](const int lIn, const int lOut, auto pF)
   {
      int dL = (lIn - lOut);
      int s = dL + pF->ls().at(0) + 2*pF->nN() - 3;

      auto p = std::make_pair(dL, s);
      return p;
   };

   // Radial function A is given
   if (pTorA && pTorB == nullptr)
   {
      this->mpOp = std::make_shared<DivR1F>(nNr + 2*q, nNc, lOut, mOut, lA, mA, lB,
         mB, pTorA, alpha, dBeta);

      auto bandInfo = prefactor(lB, lOut, pTorA);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   // Radial function B is given
   else if (pTorB && pTorA == nullptr)
   {
      this->mpOp = std::make_shared<DivR1F>(nNr + 2*q, nNc, lOut, mOut, lB, mB, lA,
         mA, pTorB, alpha, dBeta);

      auto bandInfo = prefactor(lA, lOut, pTorB);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }

   this->mIsZero = (this->elsasser(this->mLa, this->mMa, this->mLb, this->mMb,
                       this->mLout, this->mMout) == 0);
}

void ProjCurlCurlTorCrossTor::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mIsZero)
   {
      mat = Matrix::Zero(this->rows(), this->cols());
   }
   else
   {
      auto&& la = this->mLa;
      auto&& lb = this->mLb;
      auto&& lg = this->mLout;

      const Internal::MHDFloat L2g = static_cast<Internal::MHDFloat>(lg * (lg + 1));

      const MHDFloat Labg =
         this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

      Internal::MHDFloat c = L2g;

      mat = c * this->mpOp->mat();

      // Apply quasi-inverse if necessary
      this->applyQI(mat, this->mLout);

      mat *= static_cast<Internal::MHDFloat>(Labg);
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
