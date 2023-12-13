/**
 * @file IGeostrophicOperator.cpp
 * @brief Source of the implementation of the base for a geostrophic projection operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <Eigen/Dense>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "IGeostrophicOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   IGeostrophicOperator::IGeostrophicOperator(const Scalar_t ugAlpha, const Scalar_t ugDBeta, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int q)
      : IEmbeddedOperator(rows, cols, alpha, dBeta), mcUgAlpha(ugAlpha), mcUgDBeta(ugDBeta)
   {
   }

   bool IGeostrophicOperator::isUgBasis() const
   {
      return this->isUgBasis(this->mcUgAlpha, this->mcUgDBeta);
   }

   bool IGeostrophicOperator::isUgBasis(const Scalar_t a, const Scalar_t b) const
   {
      bool isBasis = (a > -1 && b > -1);

      return isBasis;
   }

   Internal::MHDFloat IGeostrophicOperator::Akj(const int k, const int j) const
   {
      Internal::MHDFloat dk = static_cast<Internal::MHDFloat>(k);
      Internal::MHDFloat dj = static_cast<Internal::MHDFloat>(j);

      Internal::MHDFloat tmp = 0;
      for(int n = k; n < j+1; n++)
      {
         Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);
         tmp += Internal::Math::pow(MHD_MP(-1),n)*Internal::Math::exp(Internal::Math::lgamma(dj+MHD_MP(1)) - Internal::Math::lgamma(dj-dn+MHD_MP(1))) / (MHD_MP(2)*dn+MHD_MP(1)) * Internal::Math::exp(Internal::Math::lgamma(MHD_MP(2)*dn+MHD_MP(2)) + Internal::Math::lgamma(dn+dk+MHD_MP(2)) - Internal::Math::lgamma(dn+MHD_MP(1)) - Internal::Math::lgamma(dn-dk+MHD_MP(1)) - Internal::Math::lgamma(MHD_MP(2)*dk+MHD_MP(2)*dn+MHD_MP(4)));
      }

      return 2*Internal::Math::PI*Internal::Math::sqrt((MHD_MP(4)*dk+MHD_MP(3)) / (MHD_MP(4)*Internal::Math::PI)) * Internal::Math::pow(MHD_MP(2),(2*k+2)) * tmp;
   }

   Internal::MHDFloat IGeostrophicOperator::Bjnab(const int j, const int n, const Internal::MHDFloat a, const Internal::MHDFloat b) const
   {
      Internal::MHDFloat dj = static_cast<Internal::MHDFloat>(j);
      Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);

      auto Cn = [&](const Internal::MHDFloat dn)
      {
         // normalisation factor
         return Internal::Math::sqrt((MHD_MP(2)*(MHD_MP(2)*dn+a+b+MHD_MP(1)))*Internal::Math::exp(Internal::Math::lgamma(dn+a+b+MHD_MP(1))+Internal::Math::lgamma(dn+MHD_MP(1))-Internal::Math::lgamma(dn+a+MHD_MP(1))-Internal::Math::lgamma(dn+b+MHD_MP(1))));
      };

      Internal::MHDFloat tmp = 0;
      for(int m = j; m < n+1; m++)
      {
         Internal::MHDFloat dm = static_cast<Internal::MHDFloat>(m);
         Internal::MHDFloat bFactor = Internal::Math::exp(Internal::Math::lgamma(dn+MHD_MP(2)) + Internal::Math::lgamma(dm+MHD_MP(2))-Internal::Math::lgamma(dn-dm+MHD_MP(1))-Internal::Math::lgamma(dm+MHD_MP(1))-Internal::Math::lgamma(dm-dj+MHD_MP(1))-Internal::Math::lgamma(dj+MHD_MP(1)));
         tmp += MHD_MP(1)/((dn+1)*(dm+1))*bFactor*Internal::Math::exp(Internal::Math::lgamma(a+b+dn+dm+MHD_MP(1)) - Internal::Math::lgamma(a+dm+MHD_MP(1))) * Internal::Math::pow(-MHD_MP(1),m-j);
      }

      Internal::MHDFloat ret = Cn(dn) * Internal::Math::exp(Internal::Math::lgamma(a+dn+MHD_MP(1)) - Internal::Math::lgamma(dn+MHD_MP(1)) - Internal::Math::lgamma(a+b+dn+MHD_MP(1))) * tmp;
      return ret;
   }

   Internal::MHDFloat IGeostrophicOperator::Bjn(const int j, const int n) const
   {
      Internal::MHDFloat a = MHD_MP(0.5);
      Internal::MHDFloat b = MHD_MP(1);
      Internal::MHDFloat dj = static_cast<Internal::MHDFloat>(j);
      Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);

      Internal::MHDFloat tmp = 0;
      for(int m = j; m < n+1; m++)
      {
         Internal::MHDFloat dm = static_cast<Internal::MHDFloat>(m);
         Internal::MHDFloat bFactor = Internal::Math::exp(Internal::Math::lgamma(dn+MHD_MP(2)) + Internal::Math::lgamma(dm+MHD_MP(2))-Internal::Math::lgamma(dn-dm+MHD_MP(1))-Internal::Math::lgamma(dm+MHD_MP(1))-Internal::Math::lgamma(dm-dj+MHD_MP(1))-Internal::Math::lgamma(dj+MHD_MP(1)));
         tmp += MHD_MP(1)/((dn+1)*(dm+1))*bFactor*Internal::Math::exp(Internal::Math::lgamma(a+b+dn+dm+MHD_MP(1)) - Internal::Math::lgamma(a+dm+MHD_MP(1))) * Internal::Math::pow(-MHD_MP(1),m-j);
      }

      Internal::MHDFloat ret = Internal::Math::sqrt((MHD_MP(2)*dn+MHD_MP(3))*(MHD_MP(4)*dn+MHD_MP(5))/(MHD_MP(8)*Internal::Math::PI*(dn+MHD_MP(1))))*Internal::Math::exp(Internal::Math::lgamma(a+dn+MHD_MP(1)) - Internal::Math::lgamma(dn+MHD_MP(1)) - Internal::Math::lgamma(a+b+dn+MHD_MP(1))) * tmp;
      return ret;
   }

   Internal::MHDFloat IGeostrophicOperator::Cnab(const int n, const Scalar_t a, const Scalar_t b) const
   {
      Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);
      Internal::MHDFloat ret = Internal::Math::sqrt(
            (MHD_MP(2)*(MHD_MP(2)*dn + a + b + MHD_MP(1)))
            )*Internal::Math::exp(
               MHD_MP(0.5)*(
               Internal::Math::lgamma(dn + a + b + MHD_MP(1))
               + Internal::Math::lgamma(dn + MHD_MP(1))
               - Internal::Math::lgamma(dn + a + MHD_MP(1))
               - Internal::Math::lgamma(dn + b + MHD_MP(1)))
               );
      return ret;
   }

   Internal::MHDFloat IGeostrophicOperator::Cn(const int n) const
   {
      Internal::MHDFloat dn = static_cast<Internal::MHDFloat>(n);
      Internal::MHDFloat ret = Internal::Math::sqrt(
            (MHD_MP(2)*dn + MHD_MP(3))*(MHD_MP(4)*dn + MHD_MP(5))/(MHD_MP(8)*Internal::Math::PI*(dn + MHD_MP(1)))
            );
      return ret;
   }

} // Worland
} // DenseSM
} // QuICC
