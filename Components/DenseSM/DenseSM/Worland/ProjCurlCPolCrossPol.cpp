/**
 * @file ProjCurlCPolCrossPol.cpp
 * @brief Source of the implementation of the projection r Curl (Curl Pol ^ Tor)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/DivR1CF.hpp"
#include "DenseSM/Worland/DivR1FC.hpp"
#include "DenseSM/Worland/ProjCurlCPolCrossPol.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Casts.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlCPolCrossPol::ProjCurlCPolCrossPol(const int nNr, const int nNc,
   const int q, const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t alpha,
   const Scalar_t dBeta) :
    IProjCrossOperator(nNr, nNc, q, lOut, mOut, lA, mA, lB, mB, alpha, dBeta)
{
   auto prefactor = [](const int lIn, const int lOut, auto pF)
   {
      int dL = (lIn - lOut);
      int s = dL + pF->ls().at(0) + 2*pF->nN() - 5;

      auto p = std::make_pair(dL, s);
      return p;
   };

   // Radial function A is given
   if (pPolA && pPolB == nullptr)
   {
      this->mpOp = std::make_shared<DivR1CF>(nNr + 2*q, nNc, lOut, mOut, lA, mA, lB,
         mB, pPolA, alpha, dBeta);

      auto bandInfo = prefactor(lB, lOut, pPolA);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   // Radial function B is given
   else if (pPolB && pPolA == nullptr)
   {
      this->mpOp = std::make_shared<DivR1FC>(nNr + 2*q, nNc, lOut, mOut, lB, mB, lA,
         mA, pPolB, alpha, dBeta);

      auto bandInfo = prefactor(lA, lOut, pPolB);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }

   this->mIsZero = (this->elsasser(this->mLa, this->mMa, this->mLb, this->mMb,
                       this->mLout, this->mMout) == 0);
}

void ProjCurlCPolCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   if (this->mIsZero)
   {
      mat = Matrix::Zero(this->rows(), this->cols());
   }
   else
   {
      auto&& lg = this->mLout;
      auto&& la = this->mLa;
      auto&& lb = this->mLb;

      const Internal::MHDFloat L2b = static_cast<Internal::MHDFloat>(lb * (lb + 1));

      const MHDFloat Labg =
         this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

      Internal::MHDFloat c = static_cast<Internal::MHDFloat>(std::pow(-1, lg + la + lb - 1)) * L2b;

      mat = c * this->mpOp->mat();

      // Apply quasi-inverse if necessary
      this->applyQI(mat, this->mLout);

      mat *= static_cast<Internal::MHDFloat>(Labg);
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
