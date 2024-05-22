/**
 * @file GeostrophicSolidBody.cpp
 * @brief Source of the implementation of the solid body with unit angular momentum operator for a geostrophic basis
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "GeostrophicSolidBody.hpp"
#include "Types/Internal/Math.hpp"
#include "Types/Internal/Literals.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/Evaluator/Set.hpp"
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "DenseSM/Worland/GeostrophicAngularMomentum.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   GeostrophicSolidBody::GeostrophicSolidBody(const int nN, const Scalar_t ugAlpha, const Scalar_t ugDBeta, const bool isGenericBasis)
      : IGeostrophicOperator(ugAlpha, ugDBeta, isGenericBasis, nN, 1)
   {
   }

   void GeostrophicSolidBody::buildOpImpl(Internal::Matrix& mat, const int rows, const int cols) const
   {
      assert(cols == 1);
      this->buildGenericOp(mat, rows);
   }

   void GeostrophicSolidBody::buildGenericOp(Internal::Matrix& mat, const int rows) const
   {
      mat = Internal::Matrix::Zero(rows, 1);

      if(this->mcIsGenericBasis)
      {
         const auto& a = this->mcUgAlpha;
         const auto& b = this->mcUgBeta;

         using namespace Internal::Literals;
         Internal::Array isg, isw;
         details::GeostrophicTools::computeGridS(isg, isw, 2*rows, a, b);
         
         Internal::Matrix ipoly;
         ipoly.resize(isg.size(), rows);
         Polynomial::Worland::Wnl wnl(a, b-1_mp);
         wnl.compute<Internal::MHDFloat>(ipoly, rows, 1, isg, isw, Polynomial::Worland::Evaluator::Set());

         mat.col(0) = (isg.transpose()*ipoly).transpose();

         // Normalize by angular momentum to get solid body with unit angular momentum
         GeostrophicAngularMomentum ag(rows, a, b, this->mcIsGenericBasis);
         auto agOp = ag.mat();

         auto agMom = (mat.col(0).transpose()*agOp.col(0)).value();
         mat.col(0) /= agMom;
      }
      else
      {
         mat.resize(1,1);
         mat(0,0) = MHD_MP(1);
      }
   }

} // Worland
} // DenseSM
} // QuICC
