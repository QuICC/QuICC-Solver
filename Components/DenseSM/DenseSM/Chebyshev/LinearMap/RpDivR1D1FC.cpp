/**
 * @file RpDivR1D1FC.cpp
 * @brief Source of the implementation of the spectral operator r^p 1/r D(f
 * (-lapl(*)))
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <cmath>

// Project includes
//
#include "DenseSM/Chebyshev/LinearMap/RpDivR1D1FC.hpp"
#include "DenseSM/Chebyshev/LinearMap/Utils/Operators.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Integrator/P.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/D.hpp"
#include "QuICC/Transform/Fft/Chebyshev/LinearMap/Projector/P.hpp"
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace DenseSM {

namespace Chebyshev {

namespace LinearMap {

RpDivR1D1FC::RpDivR1D1FC(const int nNr, const int nNc, const int p, const int lOut,
   const int mOut, const int lF, const int mF, const int lIn, const int mIn,
   std::shared_ptr<RadialTorPolFunction> pF, const Scalar_t lower,
   const Scalar_t upper) :
    ITripleHarmonicOperator(nNr, nNc, p, lOut, mOut, lF, mF, lIn, mIn, pF, lower,
       upper)
{
   if(this->mP != 4)
   {
      throw std::logic_error("Radial prefactor needs to be r^4");
   }
}

void RpDivR1D1FC::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   namespace cheb = Transform::Fft::Chebyshev::LinearMap;

   int rN =
      2 * (std::max(this->rows(), this->cols()) + this->mpF->nN() + 4 + 2);

   // Compute grid
   Internal::Array igrid, iweights;
   this->computeQuadrature(igrid, iweights, rN);
   const Internal::MHDFloat& lb = this->mcLower;
   const Internal::MHDFloat& ub = this->mcUpper;

   Matrix tA = Matrix::Identity(rN, this->cols());

   Matrix tB = Utils::evaluate(tA, this->cols(), lb, ub);
   Matrix td1B = Utils::evaluateD<1>(tA, this->cols(), lb, ub);
   Matrix td2B = Utils::evaluateD<2>(tA, this->cols(), lb, ub);
   Matrix td3B = Utils::evaluateD<3>(tA, this->cols(), lb, ub);

   Matrix f = this->mpF->evaluate(igrid, this->mLf, this->mMf).cast<MHDFloat>();
   Matrix sf = Utils::computeExpansion(f, this->mpF->nN(), lb, ub);
   Matrix d1f = this->mpF->evaluateDiff(1, igrid, this->mLf, this->mMf, lb, ub);

   const int l = this->mLin;
   const Internal::Array& r = igrid;

   tA = (r.cast<MHDFloat>().array() * d1f.array()).matrix().asDiagonal() *
           ((-r).cast<MHDFloat>().asDiagonal() *
                 (2.0 * td1B + r.cast<MHDFloat>().asDiagonal() * td2B) +
              l * (1 + l) * tB) +
        f.cast<MHDFloat>().asDiagonal() *
           (-2 * l * (1 + l) * tB +
              r.cast<MHDFloat>().asDiagonal() *
                 (td1B * (2.0 + l + l * l) -
                    r.cast<MHDFloat>().asDiagonal() *
                       (2.0 * td2B + r.cast<MHDFloat>().asDiagonal() * td3B)));

   tB = Utils::computeExpansion(tA, this->rows(), lb, ub);

   mat = tB.topRows(this->rows());
}

} // namespace LinearMap
} // namespace Chebyshev
} // namespace DenseSM
} // namespace QuICC
