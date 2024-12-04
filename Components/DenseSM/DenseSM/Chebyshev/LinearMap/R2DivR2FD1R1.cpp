/**
 * @file R2DivR2FD1R1.cpp
 * @brief Source of the implementation of the spectral operator r^2 1/r^2 f D(r *)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R2DivR2FD1R1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R2DivR2FD1R1::R2DivR2FD1R1(const int nNr, const int nNc, const int lOut, const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower, const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, lower, upper)
{}

void R2DivR2FD1R1::buildOpImpl(Internal::Matrix& mat, const int rows,
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
   cheb::Projector::D1Y1 TBwd;
   TBwd.init(sBwd);

   Matrix tmpA = Matrix::Identity(rN, this->cols());
   Matrix tmpB = Matrix::Zero(rN, this->cols());
   TBwd.transform(tmpB, tmpA);

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);

   tmpB = f.asDiagonal() * tmpB;

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
