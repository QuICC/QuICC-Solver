/**
 * @file R4DivR1F.cpp
 * @brief Source of the implementation of the spectral operator r^4 f/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/R4DivR1F.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

#include <iostream>
namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

R4DivR1F::R4DivR1F(const int nNr, const int nNc, const int lOut, const int lF, const int lIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower, const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, lF, lIn, pF, lower, upper)
{}

void R4DivR1F::buildOpImpl(Internal::Matrix& mat, const int rows,
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

   Matrix id = Matrix::Identity(this->cols(),this->cols());
   Matrix tmp = Matrix::Zero(rN,this->cols());
   TBwd.transform(tmp, id);

   auto f = this->mpF->evaluate(igrid, this->mLf);
   f = igrid.array().pow(3).matrix().asDiagonal()*f;

   tmp = f.asDiagonal() * tmp;

   id.resize(rN, this->cols());
   TFwd.transform(id, tmp);

   mat = id.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
