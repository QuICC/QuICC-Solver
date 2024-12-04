/**
 * @file R4DivR1D1DivR1D1R1F.cpp
 * @brief Source of the implementation of the spectral operator r^4 1/r D(1/r D(r f *))
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R4DivR1D1DivR1D1R1F.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R4DivR1D1DivR1D1R1F::R4DivR1D1DivR1D1R1F(const int nNr, const int nNc, const int lOut, const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower, const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, lower, upper)
{}

void R4DivR1D1DivR1D1R1F::buildOpImpl(Internal::Matrix& mat, const int rows,
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

   auto sFFwdA = std::make_shared<SetupType>(rN, 1, this->mpF->nN() + 2, pId);
   sFFwdA->setBounds(this->mcLower, this->mcUpper);
   sFFwdA->lock();
   cheb::Integrator::P TFFwdA;
   TFFwdA.init(sFFwdA);

   Matrix fA = this->mpF->evaluate(igrid, this->mLf, this->mMf);
   Matrix fB = Matrix::Zero(rN, 1);

   TFFwdA.transform(fB, fA);

   auto sFBwdA = std::make_shared<SetupType>(rN, 1, this->mpF->nN() + 2, pId);
   sFBwdA->setBounds(this->mcLower, this->mcUpper);
   sFBwdA->lock();
   cheb::Projector::D1Y1 TFBwdA;
   TFBwdA.init(sFBwdA);

   TFBwdA.transform(fA, fB);

   tmpB = fA.asDiagonal() * tmpB;

   auto sFFwdB = std::make_shared<SetupType>(rN, this->cols(), this->cols() + this->mpF->nN() + 2, pId);
   sFFwdB->setBounds(this->mcLower, this->mcUpper);
   sFFwdB->lock();
   cheb::Integrator::P TFFwdB;
   TFFwdB.init(sFFwdB);

   TFFwdB.transform(tmpA, tmpB);

   auto sFBwdB = std::make_shared<SetupType>(rN, this->cols(), this->cols() + this->mpF->nN() + 2, pId);
   sFBwdB->setBounds(this->mcLower, this->mcUpper);
   sFBwdB->lock();
   cheb::Projector::D<1> TFBwdB;
   TFBwdB.init(sFBwdB);

   Matrix tmpC = Matrix::Zero(rN, this->cols());
   TFBwdB.transform(tmpC, tmpA);

   tmpA = igrid.array().pow(2).matrix().asDiagonal()*tmpC - igrid.asDiagonal()*tmpB;

   auto sFwd = std::make_shared<SetupType>(rN, this->cols(), this->rows(), pId);
   sFwd->setBounds(this->mcLower, this->mcUpper);
   sFwd->lock();
   cheb::Integrator::P TFwd;
   TFwd.init(sFwd);

   TFwd.transform(tmpB, tmpA);

   mat = tmpB.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
