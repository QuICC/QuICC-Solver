/**
 * @file WorlandBase.cpp
 * @brief Source of the implementation of the base for the Jones-Worland polynomial base
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace Polynomial {

namespace Worland {

   const Internal::MHDFloat WorlandBase::ALPHA_CHEBYSHEV = -MHD_MP(0.5);
   const Internal::MHDFloat WorlandBase::DBETA_CHEBYSHEV = -MHD_MP(0.5);

   const Internal::MHDFloat WorlandBase::ALPHA_LEGENDRE = MHD_MP(0.0);;
   const Internal::MHDFloat WorlandBase::DBETA_LEGENDRE = -MHD_MP(0.5);

   const Internal::MHDFloat WorlandBase::ALPHA_CYLENERGY = MHD_MP(0.0);
   const Internal::MHDFloat WorlandBase::DBETA_CYLENERGY = MHD_MP(0.0);

   const Internal::MHDFloat WorlandBase::ALPHA_SPHENERGY = MHD_MP(0.0);
   const Internal::MHDFloat WorlandBase::DBETA_SPHENERGY = MHD_MP(0.5);

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

   WorlandBase::WorlandBase(const Internal::MHDFloat alpha, const Internal::MHDFloat dBeta)
      : mAlpha(alpha), mDBeta(dBeta)
   {
   }

   WorlandBase::~WorlandBase()
   {
   }

   Internal::MHDFloat WorlandBase::alpha(const int)
   {
      return this->mAlpha;
   }

   Internal::MHDFloat WorlandBase::dBeta()
   {
      return this->mDBeta;
   }

   Internal::MHDFloat WorlandBase::beta(const int l)
   {
      return Internal::MHDFloat(l) + this->dBeta();
   }

   void WorlandBase::computeW0l(Eigen::Ref<Internal::Matrix> iw0l, const int l, const Internal::MHDFloat alpha, const Internal::MHDFloat beta, const Internal::Array& igrid, ThreeTermRecurrence::NormalizerAB norm)
   {
      Internal::Array cs = norm(alpha, beta);

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

   Internal::Array WorlandBase::unitWPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(4);

      cs(0) = -(MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))*Internal::Math::sqrt((n - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/(n + a + b));
      if (n > MHD_MP(2.0))
      {
         cs(0) *= Internal::Math::sqrt((n + a + b - MHD_MP(1.0))/(MHD_MP(2.0)*n + a + b - MHD_MP(3.0)));
      }
      cs(1) = ((MHD_MP(2.0)*n + a + b)/MHD_MP(2.0))*Internal::Math::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b));
      cs(2) = (a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)))*Internal::Math::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b));
      cs(3) = Internal::Math::sqrt((MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/(n*(n + a)*(n + b)));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::unitWP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = Internal::Math::sqrt((a + b + MHD_MP(3.0))/(MHD_MP(4.0)*(a + MHD_MP(1.0))*(b + MHD_MP(1.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::unitWP0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(1);

      cs(0) = Internal::Math::sqrt(MHD_MP(2.0))*Internal::Math::exp(MHD_MP(0.5)*(Internal::Math::lgamma(a + b + MHD_MP(2.0)) - Internal::Math::lgamma(a + MHD_MP(1.0)) - Internal::Math::lgamma(b + MHD_MP(1.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   //
   // Unit Worland first derivative normalizers
   //
   Internal::Array WorlandBase::unitWDPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(4);

      cs(0) = -((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)))*Internal::Math::sqrt((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(n + a + b - MHD_MP(1.0))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))));
      cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)*n))*Internal::Math::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b - MHD_MP(1.0)));
      cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*Internal::Math::sqrt((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n + a + b - MHD_MP(1.0)));
      cs(3) = Internal::Math::sqrt((n + MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/((n + a)*(n + b)));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::unitWDP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);

      cs(2) = Internal::Math::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/(MHD_MP(2.0)*(a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b)));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::unitWDP0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(1);

      cs(0) = MHD_MP(2.0)*Internal::Math::sqrt(MHD_MP(2.0)*(a+b))*Internal::Math::exp(MHD_MP(0.5)*(Internal::Math::lgamma(a + b + MHD_MP(2.0)) - Internal::Math::lgamma(a + MHD_MP(1.0)) - Internal::Math::lgamma(b + MHD_MP(1.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   //
   // Unit Worland second derivative normalizers
   //
   Internal::Array WorlandBase::unitWD2Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(4);

      cs(0) = -((MHD_MP(2.0)*n + a + b)*(n + a + b - MHD_MP(1.0))/((MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*Internal::Math::sqrt((n+1)*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/((n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))));
      cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)))*Internal::Math::sqrt(MHD_MP(2.0)*n + a + b - MHD_MP(1.0));
      cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*Internal::Math::sqrt(MHD_MP(2.0)*n + a + b - MHD_MP(1.0));
      cs(3) = Internal::Math::sqrt((n + MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))/(n*n*(n + a)*(n + b)*(n + a + b - MHD_MP(2.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::unitWD2P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (Internal::Math::sqrt(MHD_MP(3.0))/MHD_MP(2.0))*Internal::Math::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/((a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b - MHD_MP(1.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::unitWD2P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(1);

      cs(0) = MHD_MP(8.0)*Internal::Math::sqrt((a + b)*(a + b - MHD_MP(1.0)))*Internal::Math::exp(MHD_MP(0.5)*(Internal::Math::lgamma(a + b + MHD_MP(2.0)) - Internal::Math::lgamma(a + MHD_MP(1.0)) - Internal::Math::lgamma(b + MHD_MP(1.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   //
   // Unit Worland third derivative normalizers
   //
   Internal::Array WorlandBase::unitWD3Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(4);

      cs(0) = -((MHD_MP(2.0)*n + a + b)*(n + a + b - MHD_MP(1.0))/((MHD_MP(2.0)*n + a + b - MHD_MP(2.0))))*Internal::Math::sqrt((n+2)*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))/((n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(1.0))));
      cs(1) = ((MHD_MP(2.0)*n + a + b)/(MHD_MP(2.0)));
      cs(2) = ((a*a - b*b)/(MHD_MP(2.0)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0))));
      cs(3) = Internal::Math::sqrt((n + MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b + MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(1.0))/(n*n*(n + a)*(n + b)*(n + a + b - MHD_MP(3.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::unitWD3P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = Internal::Math::sqrt((a + b + MHD_MP(1.0))*(a + b + MHD_MP(3.0))/((a + MHD_MP(1.0))*(b + MHD_MP(1.0))*(a + b - MHD_MP(2.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::unitWD3P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(1);

      cs(0) = MHD_MP(16.0)*Internal::Math::sqrt(MHD_MP(3.0)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0)))*Internal::Math::exp(MHD_MP(0.5)*(Internal::Math::lgamma(a + b + MHD_MP(2.0)) - Internal::Math::lgamma(a + MHD_MP(1.0)) - Internal::Math::lgamma(b + MHD_MP(1.0))));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural polynomial normalizer
   //
   Internal::Array WorlandBase::naturalWPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(4);

      cs(0) = -((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n*(n + a + b));
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::naturalWP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = MHD_MP(0.5);

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::naturalWP0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural first derivative normalizer
   //
   Internal::Array WorlandBase::naturalWDPnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(1.0));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::naturalWDP1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::naturalWDP0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(1);

      cs(0) = MHD_MP(2.0)*(a + b);

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural second derivative normalizer
   //
   Internal::Array WorlandBase::naturalWD2Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(2.0));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::naturalWD2P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(1.0)));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::naturalWD2P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(1);

      cs(0) = MHD_MP(4.0)*(a + b)*(a + b - MHD_MP(1.0));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural third derivative normalizer
   //
   Internal::Array WorlandBase::naturalWD3Pnab(const Internal::MHDFloat n, const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(3.0));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::naturalWD3P1ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(2.0)));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

   Internal::Array WorlandBase::naturalWD3P0ab(const Internal::MHDFloat a, const Internal::MHDFloat b)
   {
      Internal::Array cs(1);

      cs(0) = MHD_MP(8.0)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0));

      assert(!Internal::Math::isnan(cs.sum()));

      return cs;
   }

}
}
}
