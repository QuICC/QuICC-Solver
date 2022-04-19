/** 
 * @file WorlandBase.cpp
 * @brief Source of the implementation of the base for the Jones-Worland polynomial base
 */

// System includes
//
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"

// Project includes
//

namespace QuICC {

namespace Polynomial {

namespace Worland {

   const internal::MHDFloat WorlandBase::ALPHA_CHEBYSHEV = -MHD_MP(0.5);
   const internal::MHDFloat WorlandBase::DBETA_CHEBYSHEV = -MHD_MP(0.5);

   const internal::MHDFloat WorlandBase::ALPHA_LEGENDRE = MHD_MP(0.0);;
   const internal::MHDFloat WorlandBase::DBETA_LEGENDRE = -MHD_MP(0.5);

   const internal::MHDFloat WorlandBase::ALPHA_CYLENERGY = MHD_MP(0.0);
   const internal::MHDFloat WorlandBase::DBETA_CYLENERGY = MHD_MP(0.0);

   const internal::MHDFloat WorlandBase::ALPHA_SPHENERGY = MHD_MP(0.0);
   const internal::MHDFloat WorlandBase::DBETA_SPHENERGY = MHD_MP(0.5);

   WorlandBase::WorlandBase()
   {
      #if defined QUICC_WORLAND_TYPE_CHEBYSHEV
         this->mAlpha = ALPHA_CHEBYSHEV;
         this->mDBeta = DBETA_CHEBYSHEV;
      #elif defined QUICC_WORLAND_TYPE_LEGENDRE
         this->mAlpha = ALPHA_LEGENDRE;
         this->mDBeta = DBETA_LEGENDRE;
      #elif defined QUICC_WORLAND_TYPE_CYLENERGY
         this->mAlpha = ALPHA_CYLENERGY;
         this->mDBeta = DBETA_CYLENERGY;
      #elif defined QUICC_WORLAND_TYPE_SPHENERGY
         this->mAlpha = ALPHA_SPHENERGY;
         this->mDBeta = DBETA_SPHENERGY;
      #else
         #error "QUICC_WORLAND_TYPE_? is not defined"
      #endif //QUICC_WORLAND_TYPE_CHEBYSHEV
   }

   WorlandBase::WorlandBase(const internal::MHDFloat alpha, const internal::MHDFloat dBeta)
      : mAlpha(alpha), mDBeta(dBeta)
   {
   }

   WorlandBase::~WorlandBase()
   {
   }

   internal::MHDFloat WorlandBase::alpha(const int)
   {
      return this->mAlpha;
   }

   internal::MHDFloat WorlandBase::dBeta()
   {
      return this->mDBeta;
   }

   internal::MHDFloat WorlandBase::beta(const int l)
   {
      return internal::MHDFloat(l) + this->dBeta();
   }

   void WorlandBase::computeW0l(Eigen::Ref<internal::Matrix> iw0l, const int l, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, ThreeTermRecurrence::NormalizerAB norm)
   {
      internal::Array cs = norm(alpha, beta);

      if(l == 0)
      {
         iw0l.setConstant(MHD_MP(1.0));
      } else
      {
         iw0l.array() = igrid.array().pow(l);
      }

      iw0l.array() *= cs(0);
   }

   //
   // General polynomial normalizer
   //
   ThreeTermRecurrence::NormalizerNAB  WorlandBase::normWPnab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWPnab;
      #else 
         return &WorlandBase::naturalWPnab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerAB  WorlandBase::normWP1ab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWP1ab;
      #else 
         return &WorlandBase::naturalWP1ab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerAB  WorlandBase::normWP0ab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWP0ab;
      #else 
         return &WorlandBase::naturalWP0ab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerNAB  WorlandBase::normWDPnab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWDPnab;
      #else 
         return &WorlandBase::naturalWDPnab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerAB  WorlandBase::normWDP1ab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWDP1ab;
      #else 
         return &WorlandBase::naturalWDP1ab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerAB  WorlandBase::normWDP0ab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWDP0ab;
      #else 
         return &WorlandBase::naturalWDP0ab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerNAB  WorlandBase::normWD2Pnab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWD2Pnab;
      #else 
         return &WorlandBase::naturalWD2Pnab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerAB  WorlandBase::normWD2P1ab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWD2P1ab;
      #else 
         return &WorlandBase::naturalWD2P1ab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerAB  WorlandBase::normWD2P0ab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWD2P0ab;
      #else 
         return &WorlandBase::naturalWD2P0ab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerNAB  WorlandBase::normWD3Pnab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWD3Pnab;
      #else 
         return &WorlandBase::naturalWD3Pnab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerAB  WorlandBase::normWD3P1ab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWD3P1ab;
      #else 
         return &WorlandBase::naturalWD3P1ab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   ThreeTermRecurrence::NormalizerAB  WorlandBase::normWD3P0ab()
   {
      #ifdef QUICC_WORLAND_NORM_UNITY
         return &WorlandBase::unitWD3P0ab;
      #else 
         return &WorlandBase::naturalWD3P0ab;
      #endif //QUICC_WORLAND_NORM_UNITY
   }

   //
   // Unit Worland polynomial normalizers
   //

   internal::Array WorlandBase::unitWPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -(MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))*precision::sqrt((n - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/(n + a + b));
      if (n > MHD_MP(2.0))
      {
         cs(0) *= precision::sqrt((n + a + b - MHD_MP(1.0))/(MHD_MP(2.0)*n + a + b - MHD_MP(3.0)));
      }
      cs(1) = ((MHD_MP(2.0)*n + a + b)/MHD_MP(2.0))*precision::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b));
      cs(2) = (a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)))*precision::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b));
      cs(3) = precision::sqrt((MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/(n*(n + a)*(n + b)));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::unitWP1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = precision::sqrt((a + b + MHD_MP(3.0))/(MHD_MP(4.0)*(a + MHD_MP(1.0))*(b + MHD_MP(1.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::unitWP0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = precision::sqrt(MHD_MP(2.0))*precision::exp(MHD_MP(0.5)*(precision::lgamma(a + b + MHD_MP(2.0)) - precision::lgamma(a + MHD_MP(1.0)) - precision::lgamma(b + MHD_MP(1.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Unit Worland first derivative normalizers
   //
   internal::Array WorlandBase::unitWDPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)))*precision::sqrt((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(n + a + b - MHD_MP(1.0))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))));
      cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n))*precision::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b - MHD_MP(1.0)));
      cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*precision::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b - MHD_MP(1.0)));
      cs(3) = precision::sqrt((n + MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/((n + a)*(n + b)));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::unitWDP1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);

      cs(2) = precision::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/(MHD_MP(2.0)*(a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b)));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::unitWDP0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(2.0)*precision::sqrt(MHD_MP(2.0)*(a+b))*precision::exp(MHD_MP(0.5)*(precision::lgamma(a + b + MHD_MP(2.0)) - precision::lgamma(a + MHD_MP(1.0)) - precision::lgamma(b + MHD_MP(1.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Unit Worland second derivative normalizers
   //
   internal::Array WorlandBase::unitWD2Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((MHD_MP(2.0)*n + a + b)*(n + a + b - MHD_MP(1.0))/((MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*precision::sqrt((n+1)*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/((n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))));
      cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)))*precision::sqrt(MHD_MP(2.0)*n + a + b - MHD_MP(1.0));
      cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*precision::sqrt(MHD_MP(2.0)*n + a + b - MHD_MP(1.0));
      cs(3) = precision::sqrt((n + MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/(n*n*(n + a)*(n + b)*(n + a + b - MHD_MP(2.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::unitWD2P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (precision::sqrt(MHD_MP(3.0))/MHD_MP(2.0))*precision::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/((a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b - MHD_MP(1.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::unitWD2P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(8.0)*precision::sqrt((a + b)*(a + b - MHD_MP(1.0)))*precision::exp(MHD_MP(0.5)*(precision::lgamma(a + b + MHD_MP(2.0)) - precision::lgamma(a + MHD_MP(1.0)) - precision::lgamma(b + MHD_MP(1.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Unit Worland third derivative normalizers
   //
   internal::Array WorlandBase::unitWD3Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((MHD_MP(2.0)*n + a + b)*(n + a + b - MHD_MP(1.0))/((MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*precision::sqrt((n+2)*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/((n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(1.0))));
      cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)));
      cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))));
      cs(3) = precision::sqrt((n + MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n*n*(n + a)*(n + b)*(n + a + b - MHD_MP(3.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::unitWD3P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = precision::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/((a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b - MHD_MP(2.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::unitWD3P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(16.0)*precision::sqrt(MHD_MP(3.0)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0)))*precision::exp(MHD_MP(0.5)*(precision::lgamma(a + b + MHD_MP(2.0)) - precision::lgamma(a + MHD_MP(1.0)) - precision::lgamma(b + MHD_MP(1.0))));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural polynomial normalizer
   //
   internal::Array WorlandBase::naturalWPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n*(n + a + b));
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::naturalWP1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = MHD_MP(0.5);

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::naturalWP0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural first derivative normalizer
   //
   internal::Array WorlandBase::naturalWDPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(1.0));

      assert(!precision::isnan(cs.sum()));
      
      return cs;
   }

   internal::Array WorlandBase::naturalWDP1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::naturalWDP0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(2.0)*(a + b);

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural second derivative normalizer
   //
   internal::Array WorlandBase::naturalWD2Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(2.0));

      assert(!precision::isnan(cs.sum()));
      
      return cs;
   }

   internal::Array WorlandBase::naturalWD2P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(1.0)));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::naturalWD2P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(4.0)*(a + b)*(a + b - MHD_MP(1.0));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural third derivative normalizer
   //
   internal::Array WorlandBase::naturalWD3Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(3.0));

      assert(!precision::isnan(cs.sum()));
      
      return cs;
   }

   internal::Array WorlandBase::naturalWD3P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(2.0)));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array WorlandBase::naturalWD3P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(8.0)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

}
}
}
