/**
 * @file ProjCurlCurlTorCrossTor.cpp
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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlCurlTorCrossTor.hpp"
#include "DenseSM/Chebyshev/LinearMap/R4DivR1F.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCurlTorCrossTor::ProjCurlCurlTorCrossTor(const int nNr, const int nNc, const int p, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB,
   std::shared_ptr<RadialTorPolFunction> pTorA, std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t lower, const Scalar_t upper) :
    IProjCrossOperator(nNr, nNc, lOut, mOut, lA, mA, lB, mB, lower, upper)
{
   // Radial function A is given
   if(pTorA && pTorB == nullptr)
   {
      if(p == 4)
      {
         this->mpOp = std::make_shared<R4DivR1F>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pTorA, lower, upper);
      }
      else
      {
         throw std::logic_error("Radial prefactor is not implemented!");
      }
   }
   // Radial function B is given
   else if(pTorB && pTorA == nullptr)
   {
      if(p == 4)
      {
         this->mpOp = std::make_shared<R4DivR1F>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pTorB, lower, upper);
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

   this->mIsZero = (this->elsasser(this->mLa, this->mMa, this->mLb, this->mMb, this->mLout, this->mMout) == 0);
}

void ProjCurlCurlTorCrossTor::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   auto&& la = this->mLa;
   auto&& lb = this->mLb;
   auto&& lg = this->mLout;

   const MHDFloat L2g = static_cast<MHDFloat>(lg*(lg + 1));

   const MHDFloat Labg = this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

   MHDFloat c = L2g*Labg;

   mat = c*this->mpOp->mat();
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
