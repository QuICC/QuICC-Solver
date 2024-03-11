/**
 * @file ValueD1Diags.cpp
 * @brief Source of the implementation of the full sphere Worland ValueD1Diags sparse operator
 */

// System includes
//

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/ValueD1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace Chebyshev {

   ValueD1Diags::ValueD1Diags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::Stencil::ValueD1Diags(alpha, MHD_MP(-0.5), l)
   {
   }

   ValueD1Diags::ACoeff_t ValueD1Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();

      ACoeff_t num = 4.0*n*(n - 1.0)*(-3.0 + l1 + 2.0*n);
      ACoeff_t den = (-1.0 + l1 + 2.0*n)*(3.0 + 4.0*(n - 2.0)*n);

      return this->normalizeDiag(n,-2)*(num/den);
   }

   ValueD1Diags::ACoeff_t ValueD1Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();

      ACoeff_t num = -4.0*n*(l1 + 2.0*n);
      ACoeff_t den = (-1.0 + 2.0*n)*(1.0 + l1 + 2.0*n);

      return this->normalizeDiag(n,-1)*(num/den);
   }

   ValueD1Diags::ACoeff_t ValueD1Diags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Ones(n.size());

      return this->normalizeDiag(n,0)*val;
   }

} // Chebyshev
} // Stencil
} // Worland
} // SparseSM
} // QuICC
