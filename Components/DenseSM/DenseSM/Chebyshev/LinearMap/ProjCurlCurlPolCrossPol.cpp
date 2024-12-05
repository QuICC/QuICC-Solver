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
#include "DenseSM/Chebyshev/LinearMap/R4DivR1D1DivR1FD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/R4DivR1D1DivR1D1R1F.hpp"
#include "DenseSM/Chebyshev/LinearMap/R4DivR3D1R1FD1R1.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCurlPolCrossPol::ProjCurlCurlPolCrossPol(const int nNr, const int nNc, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB,
   std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t lower, const Scalar_t upper) :
    IProjCrossOperator(nNr, nNc, lOut, mOut, lA, mA, lB, mB, lower, upper)
{
   // Radial function A is given
   if(pPolA && pPolB == nullptr)
   {
      this->mpOpA = std::make_shared<R4DivR1D1DivR1FD1R1>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, lower, upper);
      this->mpOpB = std::make_shared<R4DivR1D1DivR1D1R1F>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, lower, upper);
      this->mpOpC = std::make_shared<R4DivR3D1R1FD1R1>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, lower, upper);
   }
   // Radial function B is given
   else if(pPolB && pPolA == nullptr)
   {
      this->mpOpA = std::make_shared<R4DivR1D1DivR1D1R1F>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, lower, upper);
      this->mpOpB = std::make_shared<R4DivR1D1DivR1FD1R1>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, lower, upper);
      this->mpOpC = std::make_shared<R4DivR3D1R1FD1R1>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, lower, upper);
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }
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
   MHDFloat cC = L2g;

   mat = cA*this->mpOpA->mat() + cB*this->mpOpB->mat() + cC*this->mpOpC->mat();
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
