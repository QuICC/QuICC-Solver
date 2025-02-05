/**
 * @file FR3D3.cpp
 * @brief Source of the implementation of the spectral operator  f r^3 D3(*)
 * 
 * Modified from FR2D2.cpp
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/FR3D3.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

FR3D3::FR3D3(const int nNr, 
             const int nNc, 
             const int lOut, 
             const int mOut, 
             const int lF, 
             const int mF, 
             const int lIn, 
             const int mIn,
             std::shared_ptr<RadialTorPolFunction> pF, 
             const Scalar_t lower, 
             const Scalar_t upper) :
         ITripleHarmonicOperator(nNr, nNc, lOut, mOut, lF, mF, lIn, mIn, pF, lower, upper)
{}

void FR3D3::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;
   typedef cheb::Integrator::P::SetupType SetupType;

   const auto pId = GridPurpose::SIMULATION;
   int rN = 2*(std::max(this->rows(), this->cols()) + this->mpF->nN() + 2 + 2 +(2-2) );

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);

   auto sBwd = std::make_shared<SetupType>(rN, this->cols(), this->cols(), pId);
   sBwd->setBounds(static_cast<MHDFloat>(this->mcLower), static_cast<MHDFloat>(this->mcUpper));
   sBwd->lock();
   //cheb::Projector::P TBwd;
   //TBwd.init(sBwd);

   Matrix tA = Matrix::Identity(rN,this->cols());
   //Matrix tB = Matrix::Zero(rN,this->cols());
   //TBwd.transform(tB, tA);

   //cheb::Projector::D<1> TD1Bwd;
   //TD1Bwd.init(sBwd);

   //Matrix td1B = Matrix::Zero(rN,this->cols());
   //TD1Bwd.transform(td1B, tA);

   cheb::Projector::D<3> TD3Bwd;
   TD3Bwd.init(sBwd);

   Matrix td3B = Matrix::Zero(rN,this->cols());
   TD3Bwd.transform(td3B, tA);

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();

   auto sFFwd = std::make_shared<SetupType>(rN, 1, this->mpF->nN(), pId);
   sFFwd->setBounds(static_cast<MHDFloat>(this->mcLower), static_cast<MHDFloat>(this->mcUpper));
   sFFwd->lock();
   cheb::Integrator::P TFFwd;
   TFFwd.init(sFFwd);

   Matrix sf = Matrix::Zero(rN, 1);
   TFFwd.transform(sf, f);

   //auto sFBwd = std::make_shared<SetupType>(rN, 1, this->mpF->nN(), pId);
   //sFBwd->setBounds(static_cast<MHDFloat>(this->mcLower), static_cast<MHDFloat>(this->mcUpper));
   //sFBwd->lock();
   //cheb::Projector::D<1> TFd1Bwd;
   //TFd1Bwd.init(sFBwd);

   //Matrix d1f = Matrix::Zero(rN, 1);
   //TFd1Bwd.transform(d1f, sf);

   const int l = this->mLin;
   const Internal::Array& r = igrid; 

   tA = ( f.array() * r.array().pow(3) ).cast<MHDFloat>().matrix().asDiagonal()*td3B;

   auto sFwd = std::make_shared<SetupType>(rN, this->cols(), this->rows(), pId);
   sFwd->setBounds(static_cast<MHDFloat>(this->mcLower), static_cast<MHDFloat>(this->mcUpper));
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
