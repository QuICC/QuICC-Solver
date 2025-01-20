/**
 * @file ProjCurlTorCrossTor.cpp
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
#include "DenseSM/Chebyshev/LinearMap/ProjCurlTorCrossTor.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

ProjCurlTorCrossTor::ProjCurlTorCrossTor(const int nNr, const int nNc,
   const int p, const int lOut, const int mOut, const int lA, const int mA,
   const int lB, const int mB, std::shared_ptr<RadialTorPolFunction> pTorA,
   std::shared_ptr<RadialTorPolFunction> pTorB, const Scalar_t lower,
   const Scalar_t upper) :
    IProjCrossOperator(nNr, nNc, lOut, mOut, lA, mA, lB, mB, lower, upper)
{
   // Radial function A is given
   if (pTorA && pTorB == nullptr)
   {
      // Nothing to do
   }
   // Radial function B is given
   else if (pTorB && pTorA == nullptr)
   {
      // Nothing to do
   }
   else
   {
      throw std::logic_error("One of the radial functions should be null");
   }

   this->mIsZero = true;
}

void ProjCurlTorCrossTor::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   mat = Matrix::Zero(this->rows(), this->cols());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
