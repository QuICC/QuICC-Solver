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
#include "include/QuICC/Transform/"
#include "include/QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "Types/Internal/BasicTypes.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

DivR1F::DivR1F(const int nNr, const int nNc, const int lOut, const int lF, const int lIn,
   std::shared_ptr<RadialTorPolFunction> pf, const Scalar_t lower, const Scalar_t upper) :
    ILinearMapOperator(nNr, nNc, lower, upper),
    mLout(lOut), mLf(lF), mLin(lIn), mpF(pf)
{}

void DivR1F::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   int size = 1;
   int blocks = 1;
   int specSize = 1;
   Transform::Chebyshev::LinearMap::Projector::DivY1::SetupType s(size, blocks, specSize, GridPurpose::SIMULATION);
   Transform::Chebyshev::LinearMap::Projector::DivY1 r_1Op;
   r_1Op.init();
   r_1Op.transform();

   auto f = this->mpF->evaluate(igrid, this->mLf);

   mat = opFwd.transpose() * f.asDiagonal() * opBwd;
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
