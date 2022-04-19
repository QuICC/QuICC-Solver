/** 
 * @file IDiags.cpp
 * @brief Source of the interface for Worland sparse operators
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//
#include <unsupported/Eigen/SpecialFunctions>

// Class include
//
#include "QuICC/SparseSM/Worland/IDiags.hpp"

// Project includes
//
#include "QuICC/Math/Constants.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   IDiags::IDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : mAlpha(alpha), mDBeta(dBeta), mL(static_cast<Scalar_t>(l))
   {
   }

   IDiags::~IDiags()
   {
   }

   IDiags::Scalar_t IDiags::alpha() const
   {
      return this->mAlpha;
   }

   IDiags::Scalar_t IDiags::beta(const Scalar_t l) const
   {
      return l + this->mDBeta;
   }

   IDiags::Scalar_t IDiags::l() const
   {
      return this->mL;
   }

   IDiags::ACoeff_t IDiags::lnorm(const ACoeff_t& n, const Scalar_t l) const
   {
      ACoeff_t norm;

      if(this->alpha() == 0.0)
      {
            norm = this->lnorm_legendre(n, l);
      } else if(this->alpha() == -0.5)
      {
            norm = this->lnorm_chebyshev(n, l);
      }

      return norm;
   }

   IDiags::ACoeff_t IDiags::norm(const ACoeff_t& n) const
   {
      ACoeff_t norm = this->lnorm(n, this->l());

      return norm.exp();
   }

   IDiags::ACoeff_t IDiags::invnorm(const ACoeff_t& n) const
   {
      ACoeff_t norm = -this->lnorm(n, this->l());

      return norm.exp();
   }

   IDiags::ACoeff_t IDiags::normalizeDiag(const ACoeff_t& n, const int k, const int p) const
   {
      if(n(0) + k < 0)
      {
         throw std::logic_error("Requested row normalization is inconsistent");
      }

      ACoeff_t norm = this->lnorm(n, this->l());

      Scalar_t dp = static_cast<Scalar_t>(p);
      norm -= this->lnorm(n+k, this->l() + dp);

      return norm.exp();
   }

   IDiags::ACoeff_t IDiags::lnorm_chebyshev(const ACoeff_t& n, const Scalar_t l) const
   {
      ACoeff_t norm;

      if(l == 0)
      {
         norm = -precision::log(MHD_MP(2.0)) + (n + MHD_MP(0.5)).lgamma() - (n + MHD_MP(1.0)).lgamma();
         if(n(0) == 0)
         {
            norm(0) = MHD_MP(0.5)*(precision::log(Precision::PI) - precision::log(MHD_MP(2.0)));
         }
      } else
      {
         norm = -(MHD_MP(2.0)*(MHD_MP(2.0)*n + l)).log() + (n + MHD_MP(0.5)).lgamma() + (n + l + MHD_MP(0.5)).lgamma() - (n + l).lgamma() - (n + MHD_MP(1.0)).lgamma();
         norm *= MHD_MP(0.5);
      }

      return norm;
   }

   IDiags::ACoeff_t IDiags::lnorm_legendre(const ACoeff_t& n, const Scalar_t l) const
   {
      Scalar_t b = this->beta(l);
      ACoeff_t norm = -MHD_MP(0.5)*(MHD_MP(2.0)*(MHD_MP(2.0)*n + b + MHD_MP(1.0))).log();

      return norm;
   }

}
}
}
