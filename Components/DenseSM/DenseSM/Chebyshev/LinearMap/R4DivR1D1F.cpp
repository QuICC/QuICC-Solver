/**
 * @file R4DivR1D1F.cpp
 * @brief Source of the implementation of the spectral operator r^4 1/r D(f *)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R4DivR1D1F.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

#include <iostream>
namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R4DivR1D1F::R4DivR1D1F(const int nNr, const int nNc, const int lOut, const int lF, const int lIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower, const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, lF, lIn, pF, lower, upper)
{}

void R4DivR1D1F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;

   const auto pId = GridPurpose::SIMULATION;
   int rN = 2*this->cols();

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);

   auto sBwd = std::make_shared<cheb::Projector::P::SetupType>(rN, this->cols(), this->cols(), pId);
   sBwd->setBounds(this->mcLower, this->mcUpper);
   sBwd->lock();
   cheb::Projector::P TBwd;
   TBwd.init(sBwd);

   auto sFwd = std::make_shared<cheb::Projector::P::SetupType>(rN, this->cols(), this->rows(), pId);
   sFwd->setBounds(this->mcLower, this->mcUpper);
   sFwd->lock();
   cheb::Integrator::P TFwd;
   TFwd.init(sFwd);

   Matrix id = Matrix::Identity(rN,this->cols());
   Matrix tmp = Matrix::Zero(rN,this->cols());
   TBwd.transform(tmp, id);

   auto f = this->mpF->evaluate(igrid, this->mLf);

   tmp = f.asDiagonal() * tmp;

   auto sFBwd = std::make_shared<cheb::Projector::P::SetupType>(rN, this->cols(), this->cols() + this->mpF->nN(), pId);
   sFBwd->setBounds(this->mcLower, this->mcUpper);
   sFBwd->lock();
   cheb::Projector::D<1> TFBwd;
   TFBwd.init(sFBwd);

   auto sFFwd = std::make_shared<cheb::Projector::P::SetupType>(rN, this->cols(), this->cols() + this->mpF->nN(), pId);
   sFFwd->setBounds(this->mcLower, this->mcUpper);
   sFFwd->lock();
   cheb::Integrator::P TFFwd;
   TFFwd.init(sFFwd);

   TFFwd.transform(id, tmp);
   TFBwd.transform(tmp, id);

   tmp = igrid.array().pow(3).matrix().asDiagonal()*tmp;

   TFwd.transform(id, tmp);

   mat = id.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
