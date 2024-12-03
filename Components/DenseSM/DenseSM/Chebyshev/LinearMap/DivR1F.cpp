/**
 * @file DivR1F.cpp
 * @brief Source of the implementation of the spectral operator f/r
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/DivR1F.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/DivY1.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

#include <iostream>
namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

DivR1F::DivR1F(const int nNr, const int nNc, const int lOut, const int lF, const int lIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower, const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, lOut, lF, lIn, pF, lower, upper)
{}

void DivR1F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   int size = 2*this->cols();
   int blocks = this->cols();
   int specSize = this->cols();
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, size);
   std::cerr << igrid.transpose() << std::endl;

   auto pSetup = std::make_shared<Transform::Fft::Chebyshev::LinearMap::Projector::DivY1::SetupType>(size, blocks, specSize, GridPurpose::SIMULATION);
   pSetup->setBounds(this->mcLower, this->mcUpper);
   pSetup->lock();

   Transform::Fft::Chebyshev::LinearMap::Projector::DivY1 r_1TBwd;
   r_1TBwd.init(pSetup);
   Transform::Fft::Chebyshev::LinearMap::Integrator::P TFwd;
   TFwd.init(pSetup);

   Matrix id = Matrix::Identity(size,this->cols());
   Matrix tmp = Matrix::Zero(size,this->cols());
   r_1TBwd.transform(tmp, id);
   std::cerr << tmp << std::endl;

   auto f = this->mpF->evaluate(igrid, this->mLf);

   tmp = f.asDiagonal() * tmp;

   TFwd.transform(id, tmp);

}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
