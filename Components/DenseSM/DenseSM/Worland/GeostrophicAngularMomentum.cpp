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
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   GeostrophicAngularMomentum::GeostrophicAngularMomentum(const int nN, const Scalar_t ugAlpha, const Scalar_t ugDBeta, const bool isGenericBasis)
      : IGeostrophicOperator(ugAlpha, ugDBeta, isGenericBasis, nN, 1)
   {
   }

   void GeostrophicAngularMomentum::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      assert(cols == 1);
      this->buildGenericOp(mat, rows);
   }

   void GeostrophicAngularMomentum::buildGenericOp(Internal::Matrix& mat, const int rows) const
   {
      mat = Internal::Matrix::Zero(rows, 1);

      if(this->mcIsGenericBasis)
      {
         const auto& a = this->mcUgAlpha;
         const auto& b = this->mcUgBeta;

         using namespace Internal::Literals;
         Internal::Array isg, isw;
         Polynomial::Quadrature::JacobiRule rule(0.5_mp, 0_mp);
         rule.computeQuadrature(isg, isw, 2*rows);
         isg = ((1_mp + isg.array())/2_mp).sqrt().matrix();
         isw.array() *= (Internal::Math::PI/Internal::Math::sqrt(2_mp)); // 4 pi / sqrt(32)
                                                                         //
         Internal::Matrix ipoly;
         ipoly.resize(isg.size(), rows);
         Polynomial::Worland::Wnl wnl(a, b-1_mp);
         wnl.compute<Internal::MHDFloat>(ipoly, rows, 1, isg, isw, Polynomial::Worland::Evaluator::Set());

         mat.col(0) = (isg.transpose()*ipoly).transpose();
      }
      else
      {
         mat(0,0) = MHD_MP(1);
      }
   }

} // Worland
} // DenseSM
} // QuICC
