/**
 * @file RpDivR2FD1R1.cpp
 * @brief Source of the implementation of the spectral operator r^p 1/r^p f D(r
 * *)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpDivR2FD1R1.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpDivR2FD1R1::RpDivR2FD1R1(const int nNr, const int nNc, const int p, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP < 2)
   {
      throw std::logic_error("Radial prefactor needs to be at least r^p");
   }
}

void RpDivR2FD1R1::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;
   typedef cheb::Integrator::P::SetupType SetupType;

   const auto pId = GridPurpose::SIMULATION;
   int rN =
      2 * (std::max(this->rows(), this->cols()) + this->mpF->nN() + 4 + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix tA = Matrix::Identity(rN, this->cols());
   Matrix tB = Utils::evaluateOp<cheb::Projector::D1Y1>(tA, this->cols(), lb, ub);

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);
   if(this->mP-2 > 0)
   {
      f = igrid.array().pow(this->mP-2).matrix().asDiagonal() * f;
   }

   tB = f.cast<MHDFloat>().asDiagonal() * tB;

   tA = Utils::computeExpansion(tB, this->rows(), lb, ub);
   mat = tA.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
