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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlCPolCrossTor.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCPolCrossTor::ProjCurlCPolCrossTor(const int nNr, const int nNc, const int p, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB,
   std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t lower, const Scalar_t upper) :
    ProjCurlTorCrossTor(nNr, nNc, p, lOut, mOut, lA, mA, lB, mB, pPolA, pTorB, lower, upper)
{
   // Radial function A is given
   if(pPolA && pTorB == nullptr)
   {
      // Nothing to do
   }
   // Radial function B is given
   else if(pTorB && pPolA == nullptr)
   {
      // Nothing to do
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
