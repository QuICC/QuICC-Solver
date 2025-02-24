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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlCurlCPolCrossPol.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1D1CF.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1D1FC.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR2CFD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR2D1R1FC.hpp"
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR1D1CF.hpp"
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR1D1FC.hpp"
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR2CFD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR2D1R1FC.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCurlCPolCrossPol::ProjCurlCurlCPolCrossPol(const int nNr, const int nNc,
   const int q, const int p, const int lOut, const int mOut, const int lA, const int mA,
   const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t lower,
   const Scalar_t upper) :
    ProjCurlCurlTorCrossPol(nNr, nNc, q, p, lOut, mOut, lA, mA, lB, mB, pPolA,
       pPolB, lower, upper)
{
   // Radial function A is given
   if (pPolA && pPolB == nullptr)
   {
      if (p == 4)
      {
         if(q == 0)
         {
            this->mpOpA = std::make_shared<RpDivR2CFD1R1>(nNr, nNc, p, lOut, mOut, lA,
               mA, lB, mB, pPolA, lower, upper);
            this->mpOpB = std::make_shared<RpDivR1D1CF>(nNr, nNc, p, lOut, mOut, lA,
               mA, lB, mB, pPolA, lower, upper);
         }
         else
         {
            // Disable general quasi-inverse calculation
            this->mQ = 0;

            this->mpOpA = std::make_shared<IqRpDivR2CFD1R1>(nNr, nNc, q, p, lOut, mOut, lA,
               mA, lB, mB, pPolA, lower, upper);
            this->mpOpB = std::make_shared<IqRpDivR1D1CF>(nNr, nNc, q, p, lOut, mOut, lA,
               mA, lB, mB, pPolA, lower, upper);
         }
      }
      else
      {
         throw std::logic_error("Radial prefactor is not implemented!");
      }
   }
   // Radial function B is given
   else if (pPolB && pPolA == nullptr)
   {
      if (p == 4)
      {
         if(q == 0)
         {
            this->mpOpA = std::make_shared<RpDivR2D1R1FC>(nNr, nNc, p, lOut, mOut, lB,
               mB, lA, mA, pPolB, lower, upper);
            this->mpOpB = std::make_shared<RpDivR1D1FC>(nNr, nNc, p, lOut, mOut, lB,
               mB, lA, mA, pPolB, lower, upper);
         }
         else
         {
            // Disable general quasi-inverse calculation
            this->mQ = 0;

            this->mpOpA = std::make_shared<IqRpDivR2D1R1FC>(nNr, nNc, q, p, lOut, mOut, lB,
               mB, lA, mA, pPolB, lower, upper);
            this->mpOpB = std::make_shared<IqRpDivR1D1FC>(nNr, nNc, q, p, lOut, mOut, lB,
               mB, lA, mA, pPolB, lower, upper);
         }
      }
      else
      {
         throw std::logic_error("Radial prefactor is not implemented!");
      }
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
