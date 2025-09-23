/**
 * @file ProjCurlPolCrossPol.cpp
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
#include "DenseSM/Utils.hpp"
#include "DenseSM/Chebyshev/LinearMap/ProjCurlPolCrossPol.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR2D1R1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR2FD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR2D1R1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/IqRpDivR2FD1R1.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlPolCrossPol::ProjCurlPolCrossPol(const int nNr, const int nNc,
   const int q, const int p, const int lOut, const int mOut, const int lA, const int mA,
   const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t lower,
   const Scalar_t upper) :
    IProjCrossOperator(nNr, nNc, q, p, lOut, mOut, lA, mA, lB, mB, lower, upper)
{
   // Radial function A is given
   if (pPolA && pPolB == nullptr)
   {
      if (p > 1)
      {
         if(this->mQ == 0)
         {
            this->mpOpA = std::make_shared<RpDivR2FD1R1>(nNr, nNc, p, lOut, mOut, lA,
               mA, lB, mB, pPolA, lower, upper);
            this->mpOpB = std::make_shared<RpDivR2D1R1F>(nNr, nNc, p, lOut, mOut, lA,
               mA, lB, mB, pPolA, lower, upper);
         }
         else
         {
            // Disable general quasi-inverse calculation
            this->mQ = 0;

            this->mpOpA = std::make_shared<IqRpDivR2FD1R1>(nNr, nNc, q, p, lOut, mOut, lA,
               mA, lB, mB, pPolA, lower, upper);
            this->mpOpB = std::make_shared<IqRpDivR2D1R1F>(nNr, nNc, q, p, lOut, mOut, lA,
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
      if (p > 1)
      {
         if(this->mQ == 0)
         {
            this->mpOpA = std::make_shared<RpDivR2D1R1F>(nNr, nNc, p, lOut, mOut, lB,
               mB, lA, mA, pPolB, lower, upper);
            this->mpOpB = std::make_shared<RpDivR2FD1R1>(nNr, nNc, p, lOut, mOut, lB,
               mB, lA, mA, pPolB, lower, upper);
         }
         else
         {
            // Disable general quasi-inverse calculation
            this->mQ = 0;

            this->mpOpA = std::make_shared<IqRpDivR2D1R1F>(nNr, nNc, q, p, lOut, mOut, lB,
               mB, lA, mA, pPolB, lower, upper);
            this->mpOpB = std::make_shared<IqRpDivR2FD1R1>(nNr, nNc, q, p, lOut, mOut, lB,
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

   this->mIsZero = (Utils::gaunt(this->mLa, this->mMa, this->mLb, this->mMb,
                       this->mLout, this->mMout) == 0);
}

void ProjCurlPolCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   using namespace Internal::Literals;

   auto&& lg = this->mLout;
   auto&& la = this->mLa;
   auto&& lb = this->mLb;

   const Internal::MHDFloat L2a = static_cast<Internal::MHDFloat>(la * (la + 1));
   const Internal::MHDFloat L2b = static_cast<Internal::MHDFloat>(lb * (lb + 1));
   const Internal::MHDFloat L2g = static_cast<Internal::MHDFloat>(lg * (lg + 1));

   const MHDFloat Kabg =
      Utils::gaunt(la, this->mMa, lb, this->mMb, lg, this->mMout);

   Internal::MHDFloat cA = L2a * (L2a - L2b - L2g) / 2.0_mp;
   Internal::MHDFloat cB = L2b * (L2a - L2b + L2g) / 2.0_mp;

   mat = cA * this->mpOpA->mpmat() + cB * this->mpOpB->mpmat();

   // Apply quasi-inverse if necessary
   this->applyQI(mat);

   mat *= static_cast<Internal::MHDFloat>(Kabg);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
