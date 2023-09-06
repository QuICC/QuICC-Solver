/** 
 * @file R1D1DivR1Diags.cpp
 * @brief Source of the implementation of the full sphere Worland R1D1DivR1Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/R1D1DivR1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace Chebyshev {

   R1D1DivR1Diags::R1D1DivR1Diags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::Stencil::R1D1DivR1Diags(alpha, MHD_MP(-0.5), l)
   {
   }

   R1D1DivR1Diags::ACoeff_t R1D1DivR1Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();

      ACoeff_t num = -2.0*n*(3.0 - 3.0*l1 - 8.0*n + 4.0*l1*n + 4.0*n*n);
      ACoeff_t den = (-1.0 + 2.0*n)*(-1.0 + l1 + 4.0*l1*n + 4.0*n*n);
      if(l1 == 1.0)
      {
         num(0) = 0.0;
         den(0) = 1.0;
      }

      return this->normalizeDiag(n,-1)*(num/den);
   }

   R1D1DivR1Diags::ACoeff_t R1D1DivR1Diags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Ones(n.size());

      return this->normalizeDiag(n,0)*val;
   }

} // Chebyshev
} // Stencil
} // Worland
} // SparseSM
} // QuICC
