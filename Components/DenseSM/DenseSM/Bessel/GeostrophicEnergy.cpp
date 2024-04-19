/**
 * @file GeostrophicEnergy.cpp
 * @brief Source of the implementation of the energy operator for the geostrophic basis
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "GeostrophicEnergy.hpp"
#include "DenseSM/Bessel/details/GeostrophicTools.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/Polynomial/Bessel/Generic.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   GeostrophicEnergy::GeostrophicEnergy(const int nN, const int nL, const Scalar_t sDNu)
      : IMatrixSMOperator(nL*nN, nN), mNn(nN), mNl(nL), mSDNu(sDNu)
   {
   }

   void GeostrophicEnergy::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      using namespace Internal::Literals;

      const auto& nN = this->mNn;
      const auto& nL = this->mNl;
      int nNug = details::GeostrophicTools::cylTruncNug(nN, nL);

      const auto nS = 2*details::GeostrophicTools::cylTruncNs(nN, nL);
      Internal::Array igridx, iweight;
      Polynomial::Quadrature::JacobiRule jRule(0.5_mp, 0.0_mp);
      jRule.computeQuadrature(igridx, iweight, nS);
      Internal::Array igrids =
         ((igridx.array() + 1.0_mp) / 2.0_mp).sqrt().matrix();
      iweight *= 4_mp*Internal::Math::PI / Internal::Math::sqrt(32.0_mp);

      Internal::Matrix ipoly;
      ipoly.resize(nS, nNug);
      const auto& sDNu = this->mSDNu;
      Polynomial::Bessel::Generic<Polynomial::Bessel::SphJnl> ugJnl(sDNu);
            ugJnl.compute<Internal::MHDFloat>(ipoly, nNug, 1, igrids, Internal::Array());

      mat = ipoly.transpose() * iweight.matrix().asDiagonal() * ipoly;
   }

} // Bessel
} // DenseSM
} // QuICC
