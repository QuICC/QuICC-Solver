/**
 * @file ProjCurlCPolCrossPol.cpp
 * @brief Source of the implementation of the projection r Curl (Curl Pol ^ Tor)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/ProjCurlCPolCrossPol.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1CF.hpp"
#include "DenseSM/Chebyshev/LinearMap/RpDivR1FC.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCPolCrossPol::ProjCurlCPolCrossPol(const int nNr, const int nNc,
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
         this->mpOp = std::make_shared<RpDivR1CF>(nNr, nNc, p, lOut, mOut, lA, mA,
            lB, mB, pPolA, lower, upper);
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
         this->mpOp = std::make_shared<RpDivR1FC>(nNr, nNc, p, lOut, mOut, lB, mB,
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

void ProjCurlCPolCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   auto&& lg = this->mLout;
   auto&& la = this->mLa;
   auto&& lb = this->mLb;

   const Internal::MHDFloat L2b = static_cast<Internal::MHDFloat>(lb * (lb + 1));

   const MHDFloat Labg =
      this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

   Internal::MHDFloat c = static_cast<Internal::MHDFloat>(std::pow(-1, lg + la + lb - 1)) * L2b;

   mat = c * this->mpOp->mpmat();

   // Apply quasi-inverse if necessary
   this->applyQI(mat);

   mat *= static_cast<Internal::MHDFloat>(Labg);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
