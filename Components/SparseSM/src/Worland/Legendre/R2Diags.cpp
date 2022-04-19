/** 
 * @file R2Diags.cpp
 * @brief Source of the implementation of the full sphere Worland R2Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Legendre/R2Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Legendre {

   R2Diags::R2Diags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::R2Diags(alpha, MHD_MP(-0.5), l)
   {
   }

   R2Diags::~R2Diags()
   {
   }

   R2Diags::ACoeff_t R2Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 2.0*n*(2.0*l1 + 2.0*n - 1.0)/((2.0*l1 + 4.0*n - 3.0)*(2.0*l1 + 4.0*n - 1.0));

      return this->normalizeDiag(n, -1)*val;
   }

   R2Diags::ACoeff_t R2Diags::d0(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      auto l2 = precision::pow(l1, 2);
      ACoeff_t val;

      val = (4.0*l2 + 8.0*l1*n + 8.0*n.pow(2) + 4.0*n - 1.0)/((2.0*l1 + 4.0*n - 1.0)*(2.0*l1 + 4.0*n + 3.0));

      return this->normalizeDiag(n, 0)*val;
   }

   R2Diags::ACoeff_t R2Diags::d1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t val;

      val = 2.0*(n + 1.0)*(2.0*l1 + 2.0*n + 1.0)/((2.0*l1 + 4.0*n + 3.0)*(2.0*l1 + 4.0*n + 5.0));

      return this->normalizeDiag(n, 1)*val;
   }

}
}
}
}
