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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlCurlCPolCrossTor.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1CF.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1FC.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCurlCPolCrossTor::ProjCurlCurlCPolCrossTor(const int nNr, const int nNc,
   const int q , const int p, const int lOut, const int mOut, const int lA, const int mA,
   const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t lower,
   const Scalar_t upper) :
    ProjCurlCurlTorCrossTor(nNr, nNc, q, p, lOut, mOut, lA, mA, lB, mB, pPolA,
       pTorB, lower, upper)
{
   // Radial function A is given
   if (pPolA && pTorB == nullptr)
   {
      if (p == 4)
      {
         this->mpOp = std::make_shared<RpDivR1CF>(nNr, nNc, p, lOut, mOut, lA, mA,
            lB, mB, pPolA, lower, upper);
      }
      else
      {
         throw std::logic_error("Radial prefactor is not implemented!");
      }
   }
   // Radial function B is given
   else if (pTorB && pPolA == nullptr)
   {
      this->mpOp = std::make_shared<RpDivR1FC>(nNr, nNc, p, lOut, mOut, lB, mB, lA,
         mA, pTorB, lower, upper);
   }
   else
   {
      if (p == 4)
      {
         throw std::logic_error("One of the radial functions should be null");
      }
      else
      {
         throw std::logic_error("Radial prefactor is not implemented!");
      }
   }
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
