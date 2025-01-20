/**
 * @file R4DivR2D1R1FC.cpp
 * @brief Source of the implementation of the spectral operator r^4 1/r D(r f)
 * (-lapl(*))/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R4DivR2D1R1FC.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R4DivR2D1R1FC::R4DivR2D1R1FC(const int nNr, const int nNc, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{}

void R4DivR2D1R1FC::buildOpImpl(Internal::Matrix& mat, const int rows,
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
   cheb::Projector::P TBwd;
   TBwd.init(sBwd);

   Matrix tA = Matrix::Identity(rN, this->cols());
   Matrix tB = Matrix::Zero(rN, this->cols());
   TBwd.transform(tB, tA);

   cheb::Projector::D<1> TD1Bwd;
   TD1Bwd.init(sBwd);

   Matrix td1B = Matrix::Zero(rN, this->cols());
   TD1Bwd.transform(td1B, tA);

   cheb::Projector::D<2> TD2Bwd;
   TD2Bwd.init(sBwd);

   Matrix td2B = Matrix::Zero(rN, this->cols());
   TD2Bwd.transform(td2B, tA);

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();

   auto sFFwd = std::make_shared<SetupType>(rN, 1, this->mpF->nN(), pId);
   sFFwd->setBounds(static_cast<MHDFloat>(this->mcLower),
      static_cast<MHDFloat>(this->mcUpper));
   sFFwd->lock();
   cheb::Integrator::P TFFwd;
   TFFwd.init(sFFwd);

   Matrix sf = Matrix::Zero(rN, 1);
   TFFwd.transform(sf, f);

   auto sFBwd = std::make_shared<SetupType>(rN, 1, this->mpF->nN(), pId);
   sFBwd->setBounds(static_cast<MHDFloat>(this->mcLower),
      static_cast<MHDFloat>(this->mcUpper));
   sFBwd->lock();
   cheb::Projector::D<1> TFd1Bwd;
   TFd1Bwd.init(sFBwd);

   Matrix d1f = Matrix::Zero(rN, 1);
   TFd1Bwd.transform(d1f, sf);

   const int l = this->mLin;
   const Internal::Array& r = igrid;


   tA = -((f.array() + d1f.array() * r.cast<MHDFloat>().array())
             .matrix()
             .asDiagonal() *
          (r.cast<MHDFloat>().asDiagonal() *
                (2.0 * td1B + r.cast<MHDFloat>().asDiagonal() * td2B) -
             l * (1.0 + l) * tB));

   auto sFwd = std::make_shared<SetupType>(rN, this->cols(), this->rows(), pId);
   sFwd->setBounds(static_cast<MHDFloat>(this->mcLower),
      static_cast<MHDFloat>(this->mcUpper));
   sFwd->lock();
   cheb::Integrator::P TFwd;
   TFwd.init(sFwd);

   Matrix tC = Matrix::Identity(rN, this->cols());
   TFwd.transform(tC, tA);

   mat = tC.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
