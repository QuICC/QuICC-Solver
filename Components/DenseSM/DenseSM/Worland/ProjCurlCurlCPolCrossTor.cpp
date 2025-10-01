/**
 * @file ProjCurlCurlCPolCrossTor.cpp
 * @brief Source of the implementation of the projection r Curl Curl (CurlPolA ^
 * TorB)
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
#include "DenseSM/Worland/ProjCurlCurlCPolCrossTor.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlCurlCPolCrossTor::ProjCurlCurlCPolCrossTor(const int nNr, const int nNc,
   const int q, const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ProjCurlCurlTorCrossTor(nNr, nNc, q, lOut, mOut, lA, mA, lB, mB, pPolA, pTorB,
       alpha, dBeta)
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
      this->mpOp = std::make_shared<DivR1CF>(nNr + 2*q, nNc, lOut, mOut, lA, mA, lB,
         mB, pPolA, alpha, dBeta);

      auto bandInfo = prefactor(lB, lOut, pPolA);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   // Radial function B is given
   else if (pTorB && pPolA == nullptr)
   {
      this->mpOp = std::make_shared<DivR1FC>(nNr + 2*q, nNc, lOut, mOut, lB, mB, lA,
         mA, pTorB, alpha, dBeta);

      auto bandInfo = prefactor(lA, lOut, pTorB);
      this->setBand(bandInfo.first, bandInfo.second);
      this->mBand.first -= 1;
      this->mBand.second -= 1;
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
