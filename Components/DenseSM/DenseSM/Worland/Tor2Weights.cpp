/**
 * @file Tor2Weights.cpp
 * @brief Source of the implementation of the projection operator from the
 * toroidal scalar to the geostrophic basis
 */

// System includes
//
#include <Eigen/Dense>
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "Tor2Weights.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

Tor2Weights::Tor2Weights(const int nN, const int nL, const Scalar_t alpha,
   const Scalar_t beta, const bool isTriangular) :
    IMatrixSMOperator(nN, nL * nN),
    mNn(nN),
    mNl(nL),
    mAlpha(alpha),
    mBeta(beta),
    mIsTriangular(isTriangular)
{}

void Tor2Weights::buildOpImpl(Internal::Matrix& mat, const int rows,
   const int cols) const
{
   const auto& nN = this->mNn;
   const auto& nL = this->mNl;
   const auto nS = details::GeostrophicTools::cylTruncNs(nN, nL, this->mIsTriangular);
   const auto& alphaB = this->mAlpha;
   const auto& betaB = this->mBeta;

   // compute Gauss-Jacobi quadrature in x for second alpha,beta pair
   Internal::Array igridx, ilambda;
   Polynomial::Quadrature::JacobiRule jRuleB(alphaB, betaB);
   jRuleB.computeQuadrature(igridx, ilambda, nS);

   // weights for computing misfit function
   mat = 1.0 / 8.0 / std::sqrt(2.0) * ilambda.cast<MHDFloat>();
}

} // namespace Worland
} // namespace DenseSM
} // namespace QuICC
