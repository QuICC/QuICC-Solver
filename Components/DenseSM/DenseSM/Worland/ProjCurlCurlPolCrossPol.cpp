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
#include "DenseSM/Worland/ProjCurlCurlPolCrossPol.hpp"
#include "DenseSM/Worland/DivR1D1DivR1FD1R1.hpp"
#include "DenseSM/Worland/DivR1D1DivR1D1R1F.hpp"
#include "DenseSM/Worland/DivR3D1R1FD1R1.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

ProjCurlCurlPolCrossPol::ProjCurlCurlPolCrossPol(const int nNr, const int nNc, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB,
   std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t alpha, const Scalar_t dBeta) :
    IProjCrossOperator(nNr, nNc, lOut, mOut, lA, mA, lB, mB, alpha, dBeta)
{
   // Radial function A is given
   if(pPolA && pPolB == nullptr)
   {
      this->mpOpA = std::make_shared<DivR1D1DivR1FD1R1>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, alpha, dBeta);
      this->mpOpB = std::make_shared<DivR1D1DivR1D1R1F>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, alpha, dBeta);
      this->mpOpC = std::make_shared<DivR3D1R1FD1R1>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, alpha, dBeta);
   }
   // Radial function B is given
   else if(pPolB && pPolA == nullptr)
   {
      this->mpOpA = std::make_shared<DivR1D1DivR1D1R1F>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, alpha, dBeta);
      this->mpOpB = std::make_shared<DivR1D1DivR1FD1R1>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, alpha, dBeta);
      this->mpOpC = std::make_shared<DivR3D1R1FD1R1>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, alpha, dBeta);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }

   this->mIsZero = (this->elsasser(this->mLa, this->mMa, this->mLb, this->mMb, this->mLout, this->mMout) == 0);
}

void ProjCurlCurlPolCrossPol::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   auto&& lg = this->mLout;
   auto&& la = this->mLa;
   auto&& lb = this->mLb;

   const MHDFloat L2a = static_cast<MHDFloat>(la*(la + 1));
   const MHDFloat L2b = static_cast<MHDFloat>(lb*(lb + 1));
   const MHDFloat L2g = static_cast<MHDFloat>(lg*(lg + 1));

   const MHDFloat Labg = this->elsasser(la, this->mMa, lb, this->mMb, lg, this->mMout);

   MHDFloat cA = -L2a*Labg;
   MHDFloat cB = -L2b*Labg;
   MHDFloat cC = L2g*Labg;

   mat = cA*this->mpOpA->mat() + cB*this->mpOpB->mat() + cC*this->mpOpC->mat();
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
