/**
 * @file ProjCurlTorCrossPol.cpp
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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlTorCrossPol.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1F.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlTorCrossPol::ProjCurlTorCrossPol(const int nNr, const int nNc,
   const int q, const int p, const int lOut, const int mOut, const int lA, const int mA,
   const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pTorA,
   std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t lower,
   const Scalar_t upper) :
    IProjCrossOperator(nNr, nNc, q, p, lOut, mOut, lA, mA, lB, mB, lower, upper)
{
   // Radial function A is given
   if (pTorA && pPolB == nullptr)
   {
      if (p > 1)
      {
         this->mpOp = std::make_shared<RpDivR1F>(nNr, nNc, p, lOut, mOut, lA, mA,
            lB, mB, pTorA, lower, upper);
      }
      else
      {
         throw std::logic_error("Radial prefactor is not implemented!");
      }
   }
   // Radial function B is given
   else if (pPolB && pTorA == nullptr)
   {
      if (p > 1)
      {
         this->mpOp = std::make_shared<RpDivR1F>(nNr, nNc, p, lOut, mOut, lB, mB,
            lA, mA, pPolB, lower, upper);
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

void ProjCurlTorCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   auto&& lg = this->mLout;
   auto&& la = this->mLa;
   auto&& lb = this->mLb;

   const Internal::MHDFloat L2b = static_cast<Internal::MHDFloat>(lb * (lb + 1));

   const MHDFloat Labg =
      this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

   Internal::MHDFloat c = std::pow(-1, la + lb + lg - 1) * L2b;

   mat = c * this->mpOp->mpmat();

   // Apply quasi-inverse if necessary
   this->applyQI(mat);

   mat *= static_cast<Internal::MHDFloat>(Labg);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
