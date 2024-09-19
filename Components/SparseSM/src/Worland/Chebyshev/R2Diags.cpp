/** 
 * @file R2Diags.cpp
 * @brief Source of the implementation of the full sphere Worland R2Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Chebyshev/R2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Chebyshev {

   R2Diags::R2Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::R2Diags(alpha, MHD_MP(-0.5), l, q)
   {
      if(q > 0)
      {
         throw std::logic_error("Truncation for q>0 is not implemented");
      }
   }

   R2Diags::ACoeff_t R2Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      if(l1 == MHD_MP(0))
      {
         val = n/(MHD_MP(2.0)*(l1 + MHD_MP(2.0)*n - MHD_MP(1.0)));
         val.head(1) = MHD_MP(1.0)/(l1 + MHD_MP(1.0));
      } else
      {
         val = n*(l1 + n - MHD_MP(1.0))/((l1 + MHD_MP(2.0)*n - MHD_MP(2.0))*(l1 + MHD_MP(2.0)*n - MHD_MP(1.0)));
      }

      return this->normalizeDiag(n, -1)*val;
   }

   R2Diags::ACoeff_t R2Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      if(l1 == MHD_MP(1))
      {
         val = (MHD_MP(1.0) + n)/(MHD_MP(2.0)*(n + MHD_MP(1.0)));
      } else
      {
         auto l2 = l1*l1;
         val = (MHD_MP(2.0)*l2 + MHD_MP(4.0)*l1*n - l1 + MHD_MP(4.0)*n.pow(2) - MHD_MP(1.0))/(MHD_MP(2.0)*(l1 + MHD_MP(2.0)*n - MHD_MP(1.0))*(l1 + MHD_MP(2.0)*n + MHD_MP(1.0)));
      }

      if(l1 == MHD_MP(0) || l1 == MHD_MP(1))
      {
         val.head(1) = (MHD_MP(2.0)*l1 + MHD_MP(1.0))/(MHD_MP(2.0)*(l1 + MHD_MP(1.0)));
      }

      return this->normalizeDiag(n, 0)*val;
   }

   R2Diags::ACoeff_t R2Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = (MHD_MP(2.0)*n + MHD_MP(1.0))*(MHD_MP(2.0)*l1 + MHD_MP(2.0)*n + MHD_MP(1.0))/(MHD_MP(4.0)*(l1 + MHD_MP(2.0)*n + MHD_MP(1.0))*(l1 + MHD_MP(2.0)*n + MHD_MP(2.0)));

      return this->normalizeDiag(n, 1)*val;
   }

}
}
}
}
