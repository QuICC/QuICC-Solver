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
#include "QuICC/SparseSM/Worland/SphEnergy/R2Diags.hpp"
#include "Types/Internal/Math.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace SphEnergy {

   R2Diags::R2Diags(const Scalar_t alpha, const int l, const int q)
      : QuICC::SparseSM::Worland::R2Diags(alpha, MHD_MP(0.5), l, q)
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

      val = 2.0*n*(2.0*l1 + 2.0*n + 1.0)/((2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 1.0));

      return this->normalizeDiag(n, -1)*val;
   }

   R2Diags::ACoeff_t R2Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = Internal::Math::pow(l1, 2);
      ACoeff_t val;

      val = (4.0*l2 + 8.0*l1*n + 8.0*l1 + 8.0*n.pow(2) + 12.0*n + 3.0)/((2.0*l1 + 4.0*n + 1.0)*(2.0*l1 + 4.0*n + 5.0));

      return this->normalizeDiag(n, 0)*val;
   }

   R2Diags::ACoeff_t R2Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 2.0*(n + 1.0)*(2.0*l1 + 2.0*n + 3.0)/((2.0*l1 + 4.0*n + 5.0)*(2.0*l1 + 4.0*n + 7.0));

      return this->normalizeDiag(n, 1)*val;
   }

}
}
}
}
