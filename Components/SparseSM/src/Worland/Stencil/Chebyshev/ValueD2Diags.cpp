/** 
 * @file ValueD2Diags.cpp
 * @brief Source of the implementation of the full sphere Worland ValueD2Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/ValueD2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace Chebyshev {

   ValueD2Diags::ValueD2Diags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::Stencil::ValueD2Diags(alpha, MHD_MP(-0.5), l)
   {
   }

   ValueD2Diags::ACoeff_t ValueD2Diags::d_2(const ACoeff_t& n) const
   {
      auto l1 = this->l();

      ACoeff_t num = 4.0*(-1.0 + n)*n*(-3.0 + l1 + 2.0*n)*(19.0 + 8.0*(-3.0 + n)*n + 2.0*l1*(-5.0 + 4.0*n));
      ACoeff_t den = (-1.0 + l1 + 2.0*n)*(3.0 + 4.0*(-2.0 + n)*n)*(3.0 + 8.0*(-1.0 + n)*n + l1*(-2.0 + 8.0*n));

      return this->normalizeDiag(n,-2)*(num/den);
   }

   ValueD2Diags::ACoeff_t ValueD2Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();

      ACoeff_t num = -4.0*n*(l1 + 2.0*n)*(7.0 + 8.0*n*n + l1*(2.0 + 8.0*n));
      ACoeff_t den = (-1.0 + 2.0*n)*(1.0 + l1 + 2.0*n)*(3.0 + 8.0*n*(1.0 + n) + l1*(6.0 + 8.0*n));

      return this->normalizeDiag(n,-1)*(num/den);
   }

   ValueD2Diags::ACoeff_t ValueD2Diags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Ones(n.size());

      return this->normalizeDiag(n,0)*val;
   }

} // Chebyshev
} // Stencil
} // Worland
} // SparseSM
} // QuICC
