/**
 * @file ProjCurlCurlPolCrossPol.cpp
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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlCurlPolCrossPol.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1D1DivR1D1R1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1D1DivR1FD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR3D1R1FD1R1.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCurlPolCrossPol::ProjCurlCurlPolCrossPol(const int nNr, const int nNc,
   const int q, const int p, const int lOut, const int mOut, const int lA, const int mA,
   const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pPolA,
   std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t lower,
   const Scalar_t upper) :
    IProjCrossOperator(nNr, nNc, q, p, lOut, mOut, lA, mA, lB, mB, lower, upper)
{
   // Radial function A is given
   if (pPolA && pPolB == nullptr)
   {
      if (p == 4)
      {
         this->mpOpA = std::make_shared<RpDivR1D1DivR1FD1R1>(nNr, nNc, p, lOut,
            mOut, lA, mA, lB, mB, pPolA, lower, upper);
         this->mpOpB = std::make_shared<RpDivR1D1DivR1D1R1F>(nNr, nNc, p, lOut,
            mOut, lA, mA, lB, mB, pPolA, lower, upper);
         this->mpOpC = std::make_shared<RpDivR3D1R1FD1R1>(nNr, nNc, p, lOut, mOut,
            lA, mA, lB, mB, pPolA, lower, upper);
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
         this->mpOpA = std::make_shared<RpDivR1D1DivR1D1R1F>(nNr, nNc, p, lOut,
            mOut, lB, mB, lA, mA, pPolB, lower, upper);
         this->mpOpB = std::make_shared<RpDivR1D1DivR1FD1R1>(nNr, nNc, p, lOut,
            mOut, lB, mB, lA, mA, pPolB, lower, upper);
         this->mpOpC = std::make_shared<RpDivR3D1R1FD1R1>(nNr, nNc, p, lOut, mOut,
            lB, mB, lA, mA, pPolB, lower, upper);
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

   this->mIsZero = (this->elsasser(this->mLa, this->mMa, this->mLb, this->mMb,
                       this->mLout, this->mMout) == 0);
}

void ProjCurlCurlPolCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   auto&& lg = this->mLout;
   auto&& la = this->mLa;
   auto&& lb = this->mLb;

   const Internal::MHDFloat L2a = static_cast<Internal::MHDFloat>(la * (la + 1));
   const Internal::MHDFloat L2b = static_cast<Internal::MHDFloat>(lb * (lb + 1));
   const Internal::MHDFloat L2g = static_cast<Internal::MHDFloat>(lg * (lg + 1));

   const MHDFloat Labg =
      this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

   Internal::MHDFloat cA = -L2a;
   Internal::MHDFloat cB = -L2b;
   Internal::MHDFloat cC = L2g;

   mat = cA * this->mpOpA->mpmat() + cB * this->mpOpB->mpmat() +
         cC * this->mpOpC->mpmat();

   // Apply quasi-inverse if necessary
   this->applyQI(mat);

   mat *= static_cast<Internal::MHDFloat>(Labg);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
