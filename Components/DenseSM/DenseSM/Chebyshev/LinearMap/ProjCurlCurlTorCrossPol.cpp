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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlCurlTorCrossPol.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1D1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR2D1R1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR2FD1R1.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCurlTorCrossPol::ProjCurlCurlTorCrossPol(const int nNr, const int nNc,
   const int q, const int p, const int lOut, const int mOut, const int lA, const int mA,
   const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pTorA,
   std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t lower,
   const Scalar_t upper) :
    IProjCrossOperator(nNr, nNc, q, p, lOut, mOut, lA, mA, lB, mB, lower, upper)
{
   // Radial function A is given
   if (pTorA && pPolB == nullptr)
   {
      if (p == 4)
      {
         this->mpOpA = std::make_shared<RpDivR2FD1R1>(nNr, nNc, p, lOut, mOut, lA,
            mA, lB, mB, pTorA, lower, upper);
         this->mpOpB = std::make_shared<RpDivR1D1F>(nNr, nNc, p, lOut, mOut, lA,
            mA, lB, mB, pTorA, lower, upper);
      }
      else
      {
         throw std::logic_error("Radial prefactor is not implemented!");
      }
   }
   // Radial function B is given
   else if (pPolB && pTorA == nullptr)
   {
      if (p == 4)
      {
         this->mpOpA = std::make_shared<RpDivR2D1R1F>(nNr, nNc, p, lOut, mOut, lB,
            mB, lA, mA, pPolB, lower, upper);
         this->mpOpB = std::make_shared<RpDivR1D1F>(nNr, nNc, p, lOut, mOut, lB,
            mB, lA, mA, pPolB, lower, upper);
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

   this->mIsZero = (this->gaunt(this->mLa, this->mMa, this->mLb, this->mMb,
                       this->mLout, this->mMout) == 0);
}

void ProjCurlCurlTorCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
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
      this->gaunt(la, this->mMa, lb, this->mMb, lg, this->mMout);

   Internal::MHDFloat cA = L2g * (L2a + L2b - L2g) / 2.0_mp;
   Internal::MHDFloat cB = -L2b * (L2b - L2a - L2g) / 2.0_mp;

   mat = cA * this->mpOpA->mpmat() + cB * this->mpOpB->mpmat();

   // Apply quasi-inverse if necessary
   this->applyQI(mat);

   mat *= static_cast<Internal::MHDFloat>(Kabg);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
