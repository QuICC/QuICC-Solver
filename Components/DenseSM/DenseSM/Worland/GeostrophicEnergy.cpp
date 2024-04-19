/**
 * @file GeostrophicEnergy.cpp
 * @brief Source of the implementation of the full sphere Worland projection operator onto best energy fit
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "GeostrophicEnergy.hpp"
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/Polynomial/Jacobi/Pnab.hpp"
#include "QuICC/Polynomial/Jacobi/Evaluator/Set.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   GeostrophicEnergy::GeostrophicEnergy(const int nN, const int nL, const Scalar_t ugAlpha, const Scalar_t ugBeta, const bool isGenericBasis, const bool isTriangular)
      : IGeostrophicOperator(ugAlpha, ugBeta, isGenericBasis, nL*nN, nN), mNn(nN), mNl(nL), mIsTriangular(isTriangular)
   {
   }

   void GeostrophicEnergy::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      using namespace Internal::Literals;

      const auto& nN = this->mNn;
      const auto& nL = this->mNl;
      int nNug;
      if(this->mIsTriangular)
      {
         nNug = details::GeostrophicTools::cylTruncNug(nL, this->mIsTriangular);
      }
      else
      {
         nNug = details::GeostrophicTools::cylTruncNugC(nN, nL);
      }

      if (this->mcIsGenericBasis)
      {
         const auto nS = details::GeostrophicTools::cylTruncNs(nN, nL, this->mIsTriangular);
         Internal::Array igridx, ilambda;
         Polynomial::Quadrature::JacobiRule jRule(0.5, 1.0);
         jRule.computeQuadrature(igridx, ilambda, nS);
         Internal::Array igrids =
            ((igridx.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();

         Internal::Matrix ipoly;
         ipoly.resize(nS, nNug);
         Polynomial::Jacobi::Pnab pnab;
         pnab.compute<Internal::MHDFloat>(ipoly, nNug, this->mcUgAlpha,
            this->mcUgBeta, igridx, Internal::Array(),
            Polynomial::Jacobi::Evaluator::Set());

         for (int j = 0; j < nNug; ++j)
         {
            Internal::MHDFloat Cj = DenseSM::Worland::details::GeostrophicTools::Cnab(j, this->mcUgAlpha, this->mcUgBeta);
            ipoly.col(j) = ipoly.col(j) * Cj;
         }
         Internal::Matrix tmp =
            ipoly.transpose() * ilambda.matrix().asDiagonal() * ipoly;
         tmp *= Internal::Math::PI / Internal::Math::pow(2.0_mp, 1.5_mp);
         mat = Internal::cast(tmp);
      }
      else
      {
         mat = Matrix::Identity(nNug, nNug);
      }
   }

} // Worland
} // DenseSM
} // QuICC
