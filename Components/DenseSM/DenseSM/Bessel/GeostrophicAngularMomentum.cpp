/**
 * @file GeostrophicAngularMomentum.cpp
 * @brief Source of the implementation of the angular momentum operator for a geostrophic basis
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "GeostrophicAngularMomentum.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/Polynomial/Bessel/Generic.hpp"
#include "QuICC/Polynomial/Bessel/SphJnl.hpp"

namespace QuICC {

namespace DenseSM {

namespace Bessel {

   GeostrophicAngularMomentum::GeostrophicAngularMomentum(const int nN, const Scalar_t sDNu)
      : IMatrixSMOperator(nN, 1), mNn(nN), mSDNu(sDNu)
   {
   }

   void GeostrophicAngularMomentum::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      assert(cols == 1);
      this->buildGenericOp(mat, rows);
   }

   void GeostrophicAngularMomentum::buildGenericOp(Internal::Matrix& mat, const int rows) const
   {
      using namespace Internal::Literals;
      mat = Internal::Matrix::Zero(rows, 1);

      Internal::Array isg, isw;
      Polynomial::Quadrature::JacobiRule rule(0.5_mp, 0_mp);
      rule.computeQuadrature(isg, isw, 2*rows);
      isg = ((1_mp + isg.array())/2_mp).sqrt().matrix();
      isw.array() *= (Internal::Math::PI/Internal::Math::sqrt(2_mp)); // 4 pi / sqrt(32)

      Internal::Matrix ipoly;
      ipoly.resize(isg.size(), rows);
      Polynomial::Bessel::Generic<Polynomial::Bessel::SphJnl> jnl(this->mSDNu);
      jnl.compute<Internal::MHDFloat>(ipoly, rows, 1, isg, isw);

      mat.col(0) = (isg.transpose()*ipoly).transpose();
   }

} // Bessel
} // DenseSM
} // QuICC
