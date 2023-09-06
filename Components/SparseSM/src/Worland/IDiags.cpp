/**
 * @file IDiags.cpp
 * @brief Source of the interface for Worland sparse operators
 */

// System includes
//
#include <cassert>
#include <stdexcept>
#include <unsupported/Eigen/SpecialFunctions>

// Project includes
//
#include "QuICC/SparseSM/Worland/IDiags.hpp"
#include "QuICC/Math/Constants.hpp"
#include "QuICC/SparseSM/Worland/Tools.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   IDiags::IDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : mQ(q), mAlpha(alpha), mDBeta(dBeta), mL(static_cast<Scalar_t>(l))
   {
   }

   void IDiags::zeroLast(IDiags::ACoeff_t& val, const int n) const
   {
      if(this->mQ > 0)
      {
         if(n > 0)
         {
            int size = val.size();
            val.bottomRows(std::min(size, n)).setZero();
         }
      }
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

   void IDiags::precomputeNorm(const int maxN, const int p)
   {
      Scalar_t dl = this->l() + static_cast<Scalar_t>(p);
      ACoeffI mi = ACoeffI::LinSpaced(maxN+1, 0, maxN);
      ACoeff_t m = mi.cast<Scalar_t>();
      this->mNorm[dl] = this->lnorm(m, dl);
   }

   IDiags::ACoeff_t IDiags::normalizeDiag(const ACoeff_t& n, const int k, const int p) const
   {
      if(n(0) + k < 0)
      {
         throw std::logic_error("Requested row normalization is inconsistent");
      }

      ACoeff_t norm;
      Scalar_t dp = static_cast<Scalar_t>(p);
      if(this->mNorm.count(this->l()) > 0 && this->mNorm.count(this->l()+dp) > 0)
      {
         // Check sizes are correct
         assert(this->mNorm.at(this->l()).size() > static_cast<int>(n.bottomRows(1)(0)));
         assert(this->mNorm.at(this->l()+dp).size() > static_cast<int>(n.bottomRows(1)(0)+k));

         auto n0 = static_cast<int>(n(0));
         norm = (this->mNorm.at(this->l()).segment(n0, n.size()) - this->mNorm.at(this->l()+dp).segment(n0+k, n.size())).exp();
      }
      else
      {
         norm = (this->lnorm(n, this->l()) - this->lnorm(n+k, this->l() + dp)).exp();
      }

      return norm;
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

   WorlandKind IDiags::type() const
   {
      return Tools::identifyBasis(this->mAlpha, this->mDBeta);
   }

}
}
}
