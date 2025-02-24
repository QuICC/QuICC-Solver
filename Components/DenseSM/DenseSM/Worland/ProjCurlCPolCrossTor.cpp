/**
 * @file ProjCurlCPolCrossTor.cpp
 * @brief Source of the implementation of the projection r Curl CurlPol ^ Tor
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/ProjCurlCPolCrossTor.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlCPolCrossTor::ProjCurlCPolCrossTor(const int nNr, const int nNc,
   const int q, const int lOut, const int mOut, const int lA, const int mA, const int lB,
   const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t alpha,
   const Scalar_t dBeta) :
    ProjCurlTorCrossTor(nNr, nNc, q, lOut, mOut, lA, mA, lB, mB, pPolA, pTorB,
       alpha, dBeta)
{
   // Radial function A is given
   if (pPolA && pTorB == nullptr)
   {
      // Nothing to do
   }
   // Radial function B is given
   else if (pTorB && pPolA == nullptr)
   {
      // Nothing to do
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
