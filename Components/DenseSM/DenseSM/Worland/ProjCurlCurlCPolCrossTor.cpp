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
   const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ProjCurlCurlTorCrossTor(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, pTorB,
       alpha, dBeta)
{
   // Radial function A is given
   if (pPolA && pTorB == nullptr)
   {
      this->mpOp = std::make_shared<DivR1CF>(nNr, nNc, lOut, mOut, lA, mA, lB,
         mB, pPolA, alpha, dBeta);
   }
   // Radial function B is given
   else if (pTorB && pPolA == nullptr)
   {
      this->mpOp = std::make_shared<DivR1FC>(nNr, nNc, lOut, mOut, lB, mB, lA,
         mA, pTorB, alpha, dBeta);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
