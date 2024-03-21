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
#include "DenseSM/Worland/details/GeostrophicTools.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   GeostrophicAngularMomentum::GeostrophicAngularMomentum(const Scalar_t ugAlpha, const Scalar_t ugDBeta, const int nR, const Scalar_t alpha, const Scalar_t dBeta, const int q)
      : IGeostrophicOperator(ugAlpha, ugDBeta, nR, 1, alpha, dBeta, q)
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

      if(this->isUgBasis(this->mcUgAlpha, this->mcUgDBeta))
      {
         const auto a = this->mcUgAlpha;
         const auto b = this->mcUgDBeta + MHD_MP(1);
         const int nr = rows + 1;
         Internal::Array angMom = Internal::Array::Zero(nr);

         Internal::MHDFloat piFactor = Internal::Math::pow(Internal::Math::PI,MHD_MP(1.5));
         for(int n = 0; n < nr; n++)
         {
            Internal::MHDFloat t = 0;
            for(int j = 0; j < n+1; j++)
            {
               Internal::MHDFloat dj = static_cast<Internal::MHDFloat>(j);
               t += piFactor*details::GeostrophicTools::Bjnab(j, n, a, b)*Internal::Math::exp(Internal::Math::lgamma(dj+MHD_MP(2))-Internal::Math::lgamma(dj+MHD_MP(3.5)));
            }
            angMom(n) = t;
         }

         mat = angMom.topRows(rows);
      }
      else
      {
         mat(0,0) = MHD_MP(1);
      }
   }

} // Worland
} // DenseSM
} // QuICC
