/**
 * @file JacobiRule.cpp
 * @brief Source of the Jacobi quadrature
 */

// System includes
//
#include <iostream>

// Project includes
//
#include "Types/Internal/Math.hpp"
#include "QuICC/Polynomial/Quadrature/JacobiRule.hpp"
#include "QuICC/Polynomial/Quadrature/PrueferAlgorithm.hpp"
#include "QuICC/Polynomial/Jacobi/Pnab.hpp"
#include "QuICC/Polynomial/Jacobi/dPnab.hpp"
#include "QuICC/Polynomial/Jacobi/Evaluator/Set.hpp"
#include "QuICC/Polynomial/Jacobi/JacobiAsymptotic.hpp"

namespace QuICC {

namespace Polynomial {

namespace Quadrature {

   JacobiRule::JacobiRule(const Internal::MHDFloat alpha, const Internal::MHDFloat beta)
      : mAlpha(alpha), mBeta(beta)
   {
   }

   Internal::MHDLong JacobiRule::estimateFirstZero(const int n)
   {
      Internal::MHDLong a = static_cast<Internal::MHDLong>(this->mAlpha);
      Internal::MHDLong b = static_cast<Internal::MHDLong>(this->mBeta);
      Internal::MHDLong dn = static_cast<Internal::MHDLong>(n);

      // Compute lower bound on first zero
      Internal::MHDLong B = (b - a)*((a + b + MHD_MP_LONG(6.0))*dn + MHD_MP_LONG(2.0)*(a + b));
      Internal::MHDLong A = (MHD_MP_LONG(2.0)*dn + a + b)*(dn*(MHD_MP_LONG(2.0)*dn + a + b) + MHD_MP_LONG(2.0)*(a + b + MHD_MP_LONG(2.0)));
      Internal::MHDLong D = dn*dn*Internal::Math::pow(dn + a + b + MHD_MP_LONG(1.0),MHD_MP_LONG(2.0)) + (a + MHD_MP_LONG(1.0))*(b + MHD_MP_LONG(1.0))*(dn*dn + (a + b + MHD_MP_LONG(4.0))*dn + MHD_MP_LONG(2.0)*(a + b));
      Internal::MHDLong lbound = (B - MHD_MP_LONG(4.0)*(dn - MHD_MP_LONG(1.0))*Internal::Math::sqrt(D))/A;

      return lbound;
   }

   void JacobiRule::computeQuadrature(Internal::Array& igrid, Internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      Internal::ArrayL ig(size);
      ig.setZero();
      Internal::ArrayL iw(size);
      iw.setZero();

      auto a = static_cast<Internal::MHDLong>(this->mAlpha);
      auto b = static_cast<Internal::MHDLong>(this->mBeta);
      auto n = static_cast<Internal::MHDLong>(size);
      auto factor = Jacobi::JacobiBase::normFact(n, a, b);

      #ifndef QUICC_MULTPRECISION
      constexpr int asySwitchLimit = 32;
      if(size <= asySwitchLimit)
      {
      #endif
         // Naive implementation, high error at high n due to clustering
         auto epsilon = std::numeric_limits<Internal::MHDFloat>::epsilon();
         static constexpr int ulp = 2;
         static constexpr int maxIter = 20;
         for(int k = 0; k < size; ++k)
         {
            // Asymptotic formula (WKB) - works only for positive x.
            Internal::MHDFloat r = MHD_MP(0.0);
            int h = size - k;
            auto half = floor(size/2);
            if (k < half)
            {
               std::swap(a, b);
               h = k + 1;
            }

            auto C = (2*h+a-.5)*Internal::Math::PI/(2*n+a+b+1);
            auto tanC_2 = Internal::Math::tan(.5*C);
            auto T = C + 1/((2*n+a+b+1)*(2*n+a+b+1))
               * ((.25-a*a)/tanC_2 - (.25-b*b)*tanC_2);
            r = Internal::Math::cos(T);

            if (k < floor(size/2))
            {
               std::swap(a, b);
               r = -r;
            }

            for(int it = 0; it < maxIter; ++it)
            {
               auto pab = JacobiRule::u(r, size);
               auto dPab = JacobiRule::du(r, size);
               auto delta = -pab / dPab;

               r += delta;

               auto tol = epsilon * Internal::Math::abs(r) * ulp;
               if(Internal::Math::abs(delta) <= tol)
               {
                  break;
               }
               if (it == maxIter -1) throw std::logic_error("did not converge");
            }

            // store
            ig(k) = r;
            iw(k) = JacobiRule::du(r, size);
         }

         // Convert derivative to weights
         igrid = ig.cast<Internal::MHDFloat>();
         iweights = (factor * (MHD_MP_LONG(1.0) - ig.array().square()).array().inverse()
            *iw.array().square().inverse()).cast<Internal::MHDFloat>();

      #ifndef QUICC_MULTPRECISION
      }
      else
      {
         std::array<MHDFloat,2> vals{};
         constexpr int bndSize = 10;

         auto epsilon = std::numeric_limits<MHDFloat>::epsilon();
         static constexpr std::uint32_t ulp = 4;
         static constexpr std::uint32_t maxIter = 20;
         for(int k = 0; k < size; ++k)
         {
            // Asymptotic formula (WKB) - works only for positive x.
            Internal::MHDFloat r = MHD_MP(0.0);
            std::uint32_t h = size - k;
            auto half = floor(size/2);
            if (k < half)
            {
               std::swap(a, b);
               h = k + 1;
            }

            // This approximation is technically valid only for a,b in [-1/2,1/2] and away from the boundary, see Hale and Towsend 2013
            /// \todo to enable higher beta a better initial guess needs
            /// to be implemented
            auto C = (2*h+a-.5)*Math::PI/(2*n+a+b+1);
            auto T = C + 1/((2*n+a+b+1)*(2*n+a+b+1))
               * ((.25-a*a)/std::tan(.5*C) - (.25-b*b)*std::tan(.5*C));
            r = std::cos(T);


            auto t = std::acos(r);
            for(std::uint32_t it = 0; it < maxIter; ++it)
            {
               if((k >= bndSize) && (k < size-bndSize))
               {
                  vals = Jacobi::JacobiAsy::evalPdPatInt(size, a, b, t);
               }
               else
               {
                  vals = Jacobi::JacobiAsy::evalPdPatBnd(size, a, b, t);
               }
               auto delta = -vals[0] / vals[1];

               t += delta;

               auto tol = epsilon * std::abs(t) * ulp;
               if(std::abs(delta) <= tol)
               {
                  break;
               }
               if (it == maxIter -1) throw std::logic_error("did not converge");
            }

            // eval once more to update derivative
            if((k >= bndSize) && (k < size-bndSize))
            {
               vals = Jacobi::JacobiAsy::evalPdPatInt(size, a, b, t);
            }
            else
            {
               vals = Jacobi::JacobiAsy::evalPdPatBnd(size, a, b, t);
            }
            r = std::cos(t);
            auto dPab = vals[1];
            if (k < half)
            {
               std::swap(a, b);
               r = -r;
               if (static_cast<std::uint32_t>(n) % 2 == 0)
               {
                  dPab = -dPab;
               }
            }

            // store
            ig(k) = r;
            iw(k) = 1.0/(dPab*dPab);
         }

         // Convert derivative to weights
         igrid = ig.cast<Internal::MHDFloat>();
         iweights = (factor*iw).cast<Internal::MHDFloat>();

      }
      #endif

      // Sort the grid and weights
      this->sortQuadrature(igrid, iweights);
   }

   Internal::MHDLong   JacobiRule::p(const Internal::MHDLong xi, const int diff)
   {
      // Safety asserts
      assert(diff >= 0);

      // Get p polynomial
      if(diff == 0)
      {
         return MHD_MP_LONG(1.0)-xi*xi;

      // Get first derivative of p polynomial
      } else if(diff == 1)
      {
         return MHD_MP_LONG(-2.0)*xi;

      // Get second derivative of p polynomial
      } else if(diff == 2)
      {
         return MHD_MP_LONG(-2.0);

      } else
      {
         return MHD_MP_LONG(0.0);
      }
   }

   Internal::MHDLong   JacobiRule::q(const Internal::MHDLong xi, const int diff)
   {
      // Safety asserts
      assert(diff >= 0);

      Internal::MHDLong a = static_cast<Internal::MHDLong>(this->mAlpha);
      Internal::MHDLong b = static_cast<Internal::MHDLong>(this->mBeta);

      // Get q polynomial
      if(diff == 0)
      {
         return (b - a - (a + b + MHD_MP_LONG(2.0))*xi);

      // Get first derivative of q polynomial
      } else if(diff == 1)
      {
         return -(a + b + MHD_MP_LONG(2.0));

      // Get second derivative of q polynomial
      } else if(diff == 2)
      {
         return MHD_MP_LONG(0.0);

      } else
      {
         return MHD_MP_LONG(0.0);
      }
   }

   Internal::MHDLong   JacobiRule::r(const int n, const int diff)
   {
      // Safety asserts
      assert(diff >= 0);

      Internal::MHDLong a = static_cast<Internal::MHDLong>(this->mAlpha);
      Internal::MHDLong b = static_cast<Internal::MHDLong>(this->mBeta);

      // Get r polynomial
      if(diff == 0)
      {
         Internal::MHDLong dn = static_cast<Internal::MHDLong>(n);

         return dn*(dn + a + b + MHD_MP_LONG(1.0));

      // Get first derivative of r polynomial
      } else if(diff == 1)
      {
         return MHD_MP_LONG(0.0);

      // Get second derivative of r polynomial
      } else if(diff == 2)
      {
         return MHD_MP_LONG(0.0);

      } else
      {
         return MHD_MP_LONG(0.0);
      }
   }

   Internal::MHDLong   JacobiRule::u(const Internal::MHDLong x, const int n)
   {
      Jacobi::Pnab jz;
      Internal::Array zgrid(1);
      Internal::Matrix ztmp(1,n+1);

      zgrid(0) = x;
      jz.compute<Internal::MHDFloat>(ztmp, n+1, this->mAlpha, this->mBeta, zgrid, Internal::Array(0), Jacobi::Evaluator::Set());

      return ztmp(0,n);
   }

   Internal::MHDLong   JacobiRule::du(const Internal::MHDLong x, const int n)
   {
      Jacobi::dPnab djz;
      Internal::Array zgrid(1);
      Internal::Matrix ztmp(1,n+1);

      zgrid(0) = x;
      djz.compute<Internal::MHDFloat>(ztmp, n, this->mAlpha, this->mBeta, zgrid, Internal::Array(0), Jacobi::Evaluator::Set());

      return ztmp(0,n);
   }

}
}
}
