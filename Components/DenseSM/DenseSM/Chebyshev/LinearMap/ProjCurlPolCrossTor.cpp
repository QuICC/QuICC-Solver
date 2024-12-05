/**
 * @file ProjCurlPolCrossTor.cpp
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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlPolCrossTor.hpp"
#include "DenseSM/Chebyshev/LinearMap/R2DivR1F.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlPolCrossTor::ProjCurlPolCrossTor(const int nNr, const int nNc, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB,
   std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t lower, const Scalar_t upper) :
    IProjCrossOperator(nNr, nNc, lOut, mOut, lA, mA, lB, mB, lower, upper)
{
   // Radial function A is given
   if(pPolA && pTorB == nullptr)
   {
      this->mpOp = std::make_shared<R2DivR1F>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, lower, upper);
   }
   // Radial function B is given
   else if(pTorB && pPolA == nullptr)
   {
      this->mpOp = std::make_shared<R2DivR1F>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pTorB, lower, upper);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }
}

void ProjCurlPolCrossTor::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   auto&& la = this->mLa;
   auto&& lb = this->mLb;
   auto&& lg = this->mLout;

   const MHDFloat L2a = static_cast<MHDFloat>(la*(la + 1));

   const MHDFloat Labg = this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

   MHDFloat c = L2a*Labg;

   mat = c*this->mpOp->mat();
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
