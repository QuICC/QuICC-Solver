/**
 * @file RpDivR1F.cpp
 * @brief Source of the implementation of the spectral operator r^p f/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpDivR1F.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpDivR1F::RpDivR1F(const int nNr, const int nNc, const int p, const int lOut, const int mOut,
   const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP < 2)
   {
      throw std::logic_error("Radial prefactor should be at least r^p");
   }
}

void RpDivR1F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;
   typedef cheb::Integrator::P::SetupType SetupType;

   const auto pId = GridPurpose::SIMULATION;
   int fN = this->mpF->nN() + this->mP - 1;
   int rN =
      2 * (fN + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);
   f = igrid.array().pow(this->mP-1).matrix().asDiagonal() * f;

   auto sFwdF = std::make_shared<SetupType>(rN, 1, fN, pId);
   sFwdF->setBounds(static_cast<MHDFloat>(this->mcLower),
      static_cast<MHDFloat>(this->mcUpper));
   sFwdF->lock();
   cheb::Integrator::P TFwdF;
   TFwdF.init(sFwdF);
   Matrix sF(rN,1);
   Matrix gF(rN, 1);
   gF = f.cast<MHDFloat>();

   TFwdF.transform(sF, gF);

   mat = Matrix::Zero(rows,cols);
   expansionProduct(mat, rows, cols, sF, fN);
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
