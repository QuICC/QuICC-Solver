/**
 * @file R4DivR1D1DivR1FD1R1.cpp
 * @brief Source of the implementation of the spectral operator r^4 1/r D(1/r
 * D(r f *))
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R4DivR1D1DivR1FD1R1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D1Y1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R4DivR1D1DivR1FD1R1::R4DivR1D1DivR1FD1R1(const int nNr, const int nNc,
   const int lOut, const int mOut, const int lF, const int mF, const int lIn,
   const int mIn, std::shared_ptr<RadialTorPolFunction> pF,
   const Scalar_t lower, const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{}

void R4DivR1D1DivR1FD1R1::buildOpImpl(Internal::Matrix& mat, const int rows,
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

   auto sBwd = std::make_shared<SetupType>(rN, this->cols(), this->cols(), pId);
   sBwd->setBounds(static_cast<MHDFloat>(this->mcLower),
      static_cast<MHDFloat>(this->mcUpper));
   sBwd->lock();
   cheb::Projector::D1Y1 TBwd;
   TBwd.init(sBwd);

   Matrix tmpA = Matrix::Identity(rN, this->cols());
   Matrix tmpB = Matrix::Zero(rN, this->cols());
   TBwd.transform(tmpB, tmpA);

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);

   tmpB = f.cast<MHDFloat>().asDiagonal() * tmpB;

   auto sFFwd = std::make_shared<SetupType>(rN, this->cols(),
      this->cols() + this->mpF->nN() + 2, pId);
   sFFwd->setBounds(static_cast<MHDFloat>(this->mcLower),
      static_cast<MHDFloat>(this->mcUpper));
   sFFwd->lock();
   cheb::Integrator::P TFFwd;
   TFFwd.init(sFFwd);

   TFFwd.transform(tmpA, tmpB);

   auto sFBwd = std::make_shared<SetupType>(rN, this->cols(),
      this->cols() + this->mpF->nN() + 2, pId);
   sFBwd->setBounds(static_cast<MHDFloat>(this->mcLower),
      static_cast<MHDFloat>(this->mcUpper));
   sFBwd->lock();
   cheb::Projector::D<1> TFBwd;
   TFBwd.init(sFBwd);

   Matrix tmpC = Matrix::Zero(rN, this->cols());
   TFBwd.transform(tmpC, tmpA);

   tmpA = igrid.array().pow(2).cast<MHDFloat>().matrix().asDiagonal() * tmpC -
          igrid.cast<MHDFloat>().asDiagonal() * tmpB;

   auto sFwd = std::make_shared<SetupType>(rN, this->cols(), this->rows(), pId);
   sFwd->setBounds(static_cast<MHDFloat>(this->mcLower),
      static_cast<MHDFloat>(this->mcUpper));
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
