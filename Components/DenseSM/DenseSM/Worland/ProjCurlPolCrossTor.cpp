/**
 * @file ProjCurlPolCrossTor.cpp
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
#include "DenseSM/Worland/DivR1F.hpp"
#include "DenseSM/Worland/ProjCurlPolCrossTor.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlPolCrossTor::ProjCurlPolCrossTor(const int nNr, const int nNc,
   const int q, const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
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
   if (pPolA && pTorB == nullptr)
   {
      this->mpOp = std::make_shared<DivR1F>(nNr + 2*q, nNc, lOut, mOut, lA, mA, lB,
         mB, pPolA, alpha, dBeta);

      auto bandInfo = prefactor(lB, lOut, pPolA);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   // Radial function B is given
   else if (pTorB && pPolA == nullptr)
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

void ProjCurlPolCrossTor::buildOpImpl(Internal::Matrix& mat, const int rows,
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

      const Internal::MHDFloat L2a = static_cast<Internal::MHDFloat>(la * (la + 1));

      const MHDFloat Labg =
         this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

      Internal::MHDFloat c = L2a;

      mat = c * this->mpOp->mat();

      // Apply quasi-inverse if necessary
      this->applyQI(mat, this->mLout);

      mat *= static_cast<Internal::MHDFloat>(Labg);
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
