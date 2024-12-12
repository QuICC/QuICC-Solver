/**
 * @file R4DivR1D1CF.cpp
 * @brief Source of the implementation of the spectral operator r^4 1/r D(-lapl(f) *)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R4DivR1D1CF.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R4DivR1D1CF::R4DivR1D1CF(const int nNr, const int nNc, const int lOut, const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower, const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, lower, upper)
{}

void R4DivR1D1CF::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;
   typedef cheb::Integrator::P::SetupType SetupType;

   const auto pId = GridPurpose::SIMULATION;
   int rN = 2*(std::max(this->rows(), this->cols()) + this->mpF->nN() + 4 + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);

   auto sBwd = std::make_shared<SetupType>(rN, this->cols(), this->cols(), pId);
   sBwd->setBounds(this->mcLower, this->mcUpper);
   sBwd->lock();
   cheb::Projector::P TBwd;
   TBwd.init(sBwd);

   Matrix tA = Matrix::Identity(rN,this->cols());
   Matrix tB = Matrix::Zero(rN,this->cols());
   TBwd.transform(tB, tA);

   cheb::Projector::D<1> TD1Bwd;
   TD1Bwd.init(sBwd);

   Matrix td1B = Matrix::Zero(rN,this->cols());
   TD1Bwd.transform(td1B, tA);

   cheb::Projector::D<2> TD2Bwd;
   TD2Bwd.init(sBwd);

   Matrix td2B = Matrix::Zero(rN,this->cols());
   TD2Bwd.transform(td2B, tA);

   cheb::Projector::D<3> TD3Bwd;
   TD3Bwd.init(sBwd);

   Matrix td3B = Matrix::Zero(rN,this->cols());
   TD3Bwd.transform(td3B, tA);

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf);

   auto sFFwd = std::make_shared<SetupType>(rN, 1, this->mpF->nN(), pId);
   sFFwd->setBounds(this->mcLower, this->mcUpper);
   sFFwd->lock();
   cheb::Integrator::P TFFwd;
   TFFwd.init(sFFwd);

   Matrix sf = Matrix::Zero(rN, 1);
   TFFwd.transform(sf, f);

   auto sFBwd = std::make_shared<SetupType>(rN, 1, this->mpF->nN(), pId);
   sFBwd->setBounds(this->mcLower, this->mcUpper);
   sFBwd->lock();
   cheb::Projector::D<1> TFd1Bwd;
   TFd1Bwd.init(sFBwd);

   Matrix d1f = Matrix::Zero(rN, 1);
   TFd1Bwd.transform(d1f, sf);

   cheb::Projector::D<2> TFd2Bwd;
   TFd2Bwd.init(sFBwd);

   Matrix d2f = Matrix::Zero(rN, 1);
   TFd2Bwd.transform(d2f, sf);

   cheb::Projector::D<3> TFd3Bwd;
   TFd3Bwd.init(sFBwd);

   Matrix d3f = Matrix::Zero(rN, 1);
   TFd3Bwd.transform(d3f, sf);

   const int l = this->mLf;
   const Internal::Array& r = igrid;

   tA = (-f).asDiagonal()*l*(1+l)*((-r).asDiagonal()*td1B + 2.0*tB)+r.asDiagonal()*(((-r).array()*(2.0*d1f.array() + r.array()*d2f.array())).matrix().asDiagonal()*td1B +(d1f.array()*(2+l+l*l) - r.array()*(2.0*d2f.array()+d3f.array()*r.array())).matrix().asDiagonal()*tB);

   auto sFwd = std::make_shared<SetupType>(rN, this->cols(), this->rows(), pId);
   sFwd->setBounds(this->mcLower, this->mcUpper);
   sFwd->lock();
   cheb::Integrator::P TFwd;
   TFwd.init(sFwd);

   Matrix tC = Matrix::Identity(rN,this->cols());
   TFwd.transform(tC, tA);

   mat = tC.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
