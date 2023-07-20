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
#include "IGeostrophicOperator.hpp"

namespace QuICC {

namespace DenseSM {

namespace Worland {

   IGeostrophicOperator::IGeostrophicOperator(const Scalar_t ugAlpha, const Scalar_t ugDBeta, const int rows, const int cols, const Scalar_t alpha, const Scalar_t dBeta, const int q)
      : IWorlandOperator(rows, cols, alpha, dBeta), mcUgAlpha(ugAlpha), mcUgDBeta(ugDBeta)
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

   internal::MHDFloat IGeostrophicOperator::Akj(const int k, const int j) const
   {
      internal::MHDFloat dk = static_cast<internal::MHDFloat>(k);
      internal::MHDFloat dj = static_cast<internal::MHDFloat>(j);

      internal::MHDFloat tmp = 0;
      for(int n = k; n < j+1; n++)
      {
         internal::MHDFloat dn = static_cast<internal::MHDFloat>(n);
         tmp += precision::pow(MHD_MP(-1),n)*precision::exp(precision::lgamma(dj+MHD_MP(1)) - precision::lgamma(dj-dn+MHD_MP(1))) / (MHD_MP(2)*dn+MHD_MP(1)) * precision::exp(precision::lgamma(MHD_MP(2)*dn+MHD_MP(2)) + precision::lgamma(dn+dk+MHD_MP(2)) - precision::lgamma(dn+MHD_MP(1)) - precision::lgamma(dn-dk+MHD_MP(1)) - precision::lgamma(MHD_MP(2)*dk+MHD_MP(2)*dn+MHD_MP(4)));
      }

      return 2*Precision::PI*precision::sqrt((MHD_MP(4)*dk+MHD_MP(3)) / (MHD_MP(4)*Precision::PI)) * precision::pow(MHD_MP(2),(2*k+2)) * tmp;
   }

   internal::MHDFloat IGeostrophicOperator::Bjnab(const int j, const int n, const internal::MHDFloat a, const internal::MHDFloat b) const
   {
      internal::MHDFloat dj = static_cast<internal::MHDFloat>(j);
      internal::MHDFloat dn = static_cast<internal::MHDFloat>(n);

      auto Cn = [&](const internal::MHDFloat dn)
      {
         // normalisation factor
         return precision::sqrt((MHD_MP(2)*(MHD_MP(2)*dn+a+b+MHD_MP(1)))*precision::exp(precision::lgamma(dn+a+b+MHD_MP(1))+precision::lgamma(dn+MHD_MP(1))-precision::lgamma(dn+a+MHD_MP(1))-precision::lgamma(dn+b+MHD_MP(1))));
      };

      internal::MHDFloat tmp = 0;
      for(int m = j; m < n+1; m++)
      {
         internal::MHDFloat dm = static_cast<internal::MHDFloat>(m);
         internal::MHDFloat bFactor = precision::exp(precision::lgamma(dn+MHD_MP(2)) + precision::lgamma(dm+MHD_MP(2))-precision::lgamma(dn-dm+MHD_MP(1))-precision::lgamma(dm+MHD_MP(1))-precision::lgamma(dm-dj+MHD_MP(1))-precision::lgamma(dj+MHD_MP(1)));
         tmp += MHD_MP(1)/((dn+1)*(dm+1))*bFactor*precision::exp(precision::lgamma(a+b+dn+dm+MHD_MP(1)) - precision::lgamma(a+dm+MHD_MP(1))) * precision::pow(-MHD_MP(1),m-j);
      }

      internal::MHDFloat ret = Cn(dn) * precision::exp(precision::lgamma(a+dn+MHD_MP(1)) - precision::lgamma(dn+MHD_MP(1)) - precision::lgamma(a+b+dn+MHD_MP(1))) * tmp;
      return ret;
   }

   internal::MHDFloat IGeostrophicOperator::Bjn(const int j, const int n) const
   {
      internal::MHDFloat a = MHD_MP(0.5);
      internal::MHDFloat b = MHD_MP(1);
      internal::MHDFloat dj = static_cast<internal::MHDFloat>(j);
      internal::MHDFloat dn = static_cast<internal::MHDFloat>(n);

      internal::MHDFloat tmp = 0;
      for(int m = j; m < n+1; m++)
      {
         internal::MHDFloat dm = static_cast<internal::MHDFloat>(m);
         internal::MHDFloat bFactor = precision::exp(precision::lgamma(dn+MHD_MP(2)) + precision::lgamma(dm+MHD_MP(2))-precision::lgamma(dn-dm+MHD_MP(1))-precision::lgamma(dm+MHD_MP(1))-precision::lgamma(dm-dj+MHD_MP(1))-precision::lgamma(dj+MHD_MP(1)));
         tmp += MHD_MP(1)/((dn+1)*(dm+1))*bFactor*precision::exp(precision::lgamma(a+b+dn+dm+MHD_MP(1)) - precision::lgamma(a+dm+MHD_MP(1))) * precision::pow(-MHD_MP(1),m-j);
      }

      internal::MHDFloat ret = precision::sqrt((MHD_MP(2)*dn+MHD_MP(3))*(MHD_MP(4)*dn+MHD_MP(5))/(MHD_MP(8)*Precision::PI*(dn+MHD_MP(1))))*precision::exp(precision::lgamma(a+dn+MHD_MP(1)) - precision::lgamma(dn+MHD_MP(1)) - precision::lgamma(a+b+dn+MHD_MP(1))) * tmp;
      return ret;
   }

   internal::MHDFloat IGeostrophicOperator::Cnab(const int n, const Scalar_t a, const Scalar_t b) const
   {
      internal::MHDFloat dn = static_cast<internal::MHDFloat>(n);
      internal::MHDFloat ret = precision::sqrt(
            (MHD_MP(2)*(MHD_MP(2)*dn + a + b + MHD_MP(1)))
            )*precision::exp(
               MHD_MP(0.5)*(
               precision::lgamma(dn + a + b + MHD_MP(1))
               + precision::lgamma(dn + MHD_MP(1))
               - precision::lgamma(dn + a + MHD_MP(1))
               - precision::lgamma(dn + b + MHD_MP(1)))
               );
      return ret;
   }

   internal::MHDFloat IGeostrophicOperator::Cn(const int n) const
   {
      internal::MHDFloat dn = static_cast<internal::MHDFloat>(n);
      internal::MHDFloat ret = precision::sqrt(
            (MHD_MP(2)*dn + MHD_MP(3))*(MHD_MP(4)*dn + MHD_MP(5))/(MHD_MP(8)*Precision::PI*(dn + MHD_MP(1)))
            );
      return ret;
   }

} // Worland
} // DenseSM
} // QuICC
