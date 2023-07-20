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

namespace QuICC {

namespace DenseSM {

namespace Worland {

   GeostrophicAngularMomentum::GeostrophicAngularMomentum(const Scalar_t ugAlpha, const Scalar_t ugDBeta, const int nR, const Scalar_t alpha, const Scalar_t dBeta, const int q)
      : IGeostrophicOperator(ugAlpha, ugDBeta, nR, 1, alpha, dBeta, q)
   {
   }

   void GeostrophicAngularMomentum::buildOpImpl(internal::Matrix& mat, const int rows, const int cols) const
   {
      assert(cols == 1);
      switch(this->type())
      {
         case WorlandKind::CHEBYSHEV:
            this->buildChebyshevOp(mat, rows);
            break;
         case WorlandKind::LEGENDRE:
            throw std::logic_error("Legendre basis operator not implemented");
            break;
         case WorlandKind::CYLENERGY:
            throw std::logic_error("Cylindrical energy basis operator not implemented");
            break;
         case WorlandKind::SPHENERGY:
            throw std::logic_error("Spherical energy basis operator not implemented");
            break;
      }
   }

   void GeostrophicAngularMomentum::buildChebyshevOp(internal::Matrix& mat, const int rows) const
   {
      mat = internal::Matrix::Zero(rows, 1);

      if(this->isUgBasis(this->mcUgAlpha, this->mcUgDBeta))
      {
         const auto a = this->mcUgAlpha;
         const auto b = this->mcUgDBeta + MHD_MP(1);
         const int nr = rows + 1;
         internal::Array angMom = internal::Array::Zero(nr);

         internal::MHDFloat piFactor = precision::pow(Precision::PI,MHD_MP(1.5));
         for(int n = 0; n < nr; n++)
         {
            internal::MHDFloat t = 0;
            for(int j = 0; j < n+1; j++)
            {
               internal::MHDFloat dj = static_cast<internal::MHDFloat>(j);
               t += piFactor*this->Bjnab(j, n, a, b)*precision::exp(precision::lgamma(dj+MHD_MP(2))-precision::lgamma(dj+MHD_MP(3.5)));
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
