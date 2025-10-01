/**
 * @file ProjCurlCurlCPolCrossPol.cpp
 * @brief Source of the implementation of the projection r Curl Curl (CurlPolA ^
 * PolB)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/DivR1D1CF.hpp"
#include "DenseSM/Worland/DivR1D1FC.hpp"
#include "DenseSM/Worland/DivR2CFD1R1.hpp"
#include "DenseSM/Worland/DivR2D1R1FC.hpp"
#include "DenseSM/Worland/ProjCurlCurlCPolCrossPol.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlCurlCPolCrossPol::ProjCurlCurlCPolCrossPol(const int nNr, const int nNc,
   const int q, const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ProjCurlCurlTorCrossPol(nNr, nNc, q, lOut, mOut, lA, mA, lB, mB, pPolA, pPolB,
       alpha, dBeta)
{

   auto prefactor = [](const int lIn, const int lOut, auto pF)
   {
      int dL = (lIn - lOut);
      int s = dL + pF->ls().at(0) + 2*pF->nN() - 6;

      auto p = std::make_pair(dL, s);
      return p;
   };

   // Radial function A is given
   if (pPolA && pPolB == nullptr)
   {
      this->mpOpA = std::make_shared<DivR2CFD1R1>(nNr + 2*q, nNc, lOut, mOut, lA, mA,
         lB, mB, pPolA, alpha, dBeta);
      this->mpOpB = std::make_shared<DivR1D1CF>(nNr + 2*q, nNc, lOut, mOut, lA, mA,
         lB, mB, pPolA, alpha, dBeta);

      // Special condition where both operators cancel each other
      if (pPolA->ls().size() == 1 && pPolA->ls().at(0) == 1 && pPolA->nN() == 2)
      {
         this->mIsZero = true;
      }

      auto bandInfo = prefactor(lB, lOut, pPolA);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   // Radial function B is given
   else if (pPolB && pPolA == nullptr)
   {
      this->mpOpA = std::make_shared<DivR2D1R1FC>(nNr + 2*q, nNc, lOut, mOut, lB, mB,
         lA, mA, pPolB, alpha, dBeta);
      this->mpOpB = std::make_shared<DivR1D1FC>(nNr + 2*q, nNc, lOut, mOut, lB, mB,
         lA, mA, pPolB, alpha, dBeta);

      auto bandInfo = prefactor(lA, lOut, pPolB);
      this->setBand(bandInfo.first, bandInfo.second);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
