/**
 * @file ProjCurlCurlCPolCrossPol.cpp
 * @brief Source of the implementation of the projection r Curl Curl (CurlPolA ^ PolB)
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
#include "DenseSM/Chebyshev/LinearMap/R4DivR2CFD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/R4DivR1D1CF.hpp"
#include "DenseSM/Chebyshev/LinearMap/R4DivR2D1R1FC.hpp"
#include "DenseSM/Chebyshev/LinearMap/R4DivR1D1FC.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlCurlCPolCrossPol::ProjCurlCurlCPolCrossPol(const int nNr, const int nNc, const int lOut, const int mOut, const int lA, const int mA, const int lB, const int mB,
   std::shared_ptr<RadialTorPolFunction> pPolA, std::shared_ptr<RadialTorPolFunction> pPolB, const Scalar_t lower, const Scalar_t upper) :
    ProjCurlCurlTorCrossPol(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, pPolB, lower, upper)
{
   // Radial function A is given
   if(pPolA && pPolB == nullptr)
   {
      this->mpOpA = std::make_shared<R4DivR2CFD1R1>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, lower, upper);
      this->mpOpB = std::make_shared<R4DivR1D1CF>(nNr, nNc, lOut, mOut, lA, mA, lB, mB, pPolA, lower, upper);
   }
   // Radial function B is given
   else if(pPolB && pPolA == nullptr)
   {
      this->mpOpA = std::make_shared<R4DivR2D1R1FC>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, lower, upper);
      this->mpOpB = std::make_shared<R4DivR1D1FC>(nNr, nNc, lOut, mOut, lB, mB, lA, mA, pPolB, lower, upper);
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
