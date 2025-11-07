/**
 * @file RpFD4.cpp
 * @brief Source of the implementation of the spectral operator r^p f D2(*)
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpFD4.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpFD4::RpFD4(const int nNr, const int nNc, const int p, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP <1)
   {
      throw std::logic_error("Radial prefactor needs to be r^p, with p>0");
   }
}

void RpFD4::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;

    int rN = 2 * (std::max(this->rows(), this->cols()) + this->mpF->nN() + 4 + 2); // same as RpDivR1FC


   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix tA = Matrix::Identity(rN, this->cols());
   Matrix td4B = Utils::evaluateD<4>(tA, this->cols(), lb, ub);

   auto f = this->mpF->evaluate(igrid, this->mLf, this->mMf);
   // constructor forces mP>0
   f = igrid.array().pow(this->mP).matrix().asDiagonal() * f;
   

   td4B = f.cast<MHDFloat>().asDiagonal() * td4B;

   tA = Utils::computeExpansion(td4B, this->rows(), lb, ub);
   mat = tA.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
