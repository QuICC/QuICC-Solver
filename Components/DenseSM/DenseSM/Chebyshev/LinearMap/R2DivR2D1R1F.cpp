/**
 * @file R2DivR2D1R1F.cpp
 * @brief Source of the implementation of the spectral operator r^2 1/r^2 D(r f) *
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R2DivR2D1R1F.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R2DivR2D1R1F::R2DivR2D1R1F(const int nNr, const int nNc, const int lOut, const int lF, const int lIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower, const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, lF, lIn, pF, lower, upper)
{}

void R2DivR2D1R1F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;
   typedef cheb::Integrator::P::SetupType SetupType;

   const auto pId = GridPurpose::SIMULATION;
   int rN = 2*this->cols();

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);

   auto sBwd = std::make_shared<SetupType>(rN, this->cols(), this->cols(), pId);
   sBwd->setBounds(this->mcLower, this->mcUpper);
   sBwd->lock();
   cheb::Projector::P TBwd;
   TBwd.init(sBwd);

   Matrix tmpA = Matrix::Identity(rN, this->cols());
   Matrix tmpB = Matrix::Zero(rN, this->cols());
   TBwd.transform(tmpB, tmpA);

   auto sFFwd = std::make_shared<SetupType>(rN, 1, this->mpF->nN(), pId);
   sFFwd->setBounds(this->mcLower, this->mcUpper);
   sFFwd->lock();
   cheb::Integrator::P TFFwd;
   TFFwd.init(sFFwd);

   Matrix fB(rN,1);
   Matrix fA = this->mpF->evaluate(igrid, this->mLf);
   TFFwd.transform(fB, fA);

   auto sFBwd = std::make_shared<SetupType>(rN, 1, this->mpF->nN(), pId);
   sFBwd->setBounds(this->mcLower, this->mcUpper);
   sFBwd->lock();
   cheb::Projector::D1Y1 TFBwd;
   TFBwd.init(sFBwd);

   TFBwd.transform(fA, fB);

   tmpB = fA.asDiagonal() * tmpB;

   auto sFwd = std::make_shared<SetupType>(rN, this->cols(), this->rows(), pId);
   sFwd->setBounds(this->mcLower, this->mcUpper);
   sFwd->lock();
   cheb::Integrator::P TFwd;
   TFwd.init(sFwd);

   TFwd.transform(tmpA, tmpB);

   mat = tmpA.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
