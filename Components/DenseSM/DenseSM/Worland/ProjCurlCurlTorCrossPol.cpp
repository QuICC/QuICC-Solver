/**
 * @file ProjCurlCurlTorCrossPol.cpp
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
#include "DenseSM/Worland/DivR1D1F.hpp"
#include "DenseSM/Worland/DivR2D1R1F.hpp"
#include "DenseSM/Worland/DivR2FD1R1.hpp"
#include "DenseSM/Worland/ProjCurlCurlTorCrossPol.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlCurlTorCrossPol::ProjCurlCurlTorCrossPol(const int nNr, const int nNc,
   const int q, const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, std::shared_ptr<RadialTorPolFunction> pTorA,
   std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t alpha,
   const Scalar_t dBeta) :
    IProjCrossOperator(nNr, nNc, q, lOut, mOut, lA, mA, lB, mB, alpha, dBeta)
{
   auto prefactor = [](const int lIn, const int lOut, auto pF)
   {
      int dL = (lIn - lOut);
      int s = dL + pF->ls().at(0) + 2*pF->nN() - 4;

      auto p = std::make_pair(dL, s);
      return p;
   };

   // Radial function A is given
   if (pTorA && pPolB == nullptr)
   {
      this->mpOpA = std::make_shared<DivR2FD1R1>(nNr + 2*q, nNc, lOut, mOut, lA, mA,
         lB, mB, pTorA, alpha, dBeta);
      this->mpOpB = std::make_shared<DivR1D1F>(nNr + 2*q, nNc, lOut, mOut, lA, mA, lB,
         mB, pTorA, alpha, dBeta);

      auto bandInfo = prefactor(lB, lOut, pTorA);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   // Radial function B is given
   else if (pPolB && pTorA == nullptr)
   {
      this->mpOpA = std::make_shared<DivR2D1R1F>(nNr + 2*q, nNc, lOut, mOut, lB, mB,
         lA, mA, pPolB, alpha, dBeta);
      this->mpOpB = std::make_shared<DivR1D1F>(nNr + 2*q, nNc, lOut, mOut, lB, mB, lA,
         mA, pPolB, alpha, dBeta);

      auto bandInfo = prefactor(lA, lOut, pPolB);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }

   this->mIsZero = (this->gaunt(this->mLa, this->mMa, this->mLb, this->mMb,
                       this->mLout, this->mMout) == 0);
}

void ProjCurlCurlTorCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mIsZero)
   {
      mat = Matrix::Zero(this->rows(), this->cols());
   }
   else
   {
      using namespace Internal::Literals;
      auto&& lg = this->mLout;
      auto&& la = this->mLa;
      auto&& lb = this->mLb;

      const Internal::MHDFloat L2a = static_cast<Internal::MHDFloat>(la * (la + 1));
      const Internal::MHDFloat L2b = static_cast<Internal::MHDFloat>(lb * (lb + 1));
      const Internal::MHDFloat L2g = static_cast<Internal::MHDFloat>(lg * (lg + 1));

      const MHDFloat Kabg =
         this->gaunt(la, this->mMa, lb, this->mMb, lg, this->mMout);

      Internal::MHDFloat cA = L2g * (L2a + L2b - L2g) / 2.0_mp;
      Internal::MHDFloat cB = -L2b * (L2b - L2a - L2g) / 2.0_mp;

      mat = cA * this->mpOpA->mat() + cB * this->mpOpB->mat();

      // Apply quasi-inverse if necessary
      this->applyQI(mat, this->mLout);

      mat *= static_cast<Internal::MHDFloat>(Kabg);
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
