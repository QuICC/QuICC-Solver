/** 
 * @file D1Diags.cpp
 * @brief Source of the implementation of the full sphere Worland D1Diags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/D1Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace Chebyshev {

   D1Diags::D1Diags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::Stencil::D1Diags(alpha, MHD_MP(-0.5), l)
   {
   }

   D1Diags::ACoeff_t D1Diags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();
      ACoeff_t num = -2.0*n*(4.0*(-1.0 + n)*(-1.0 + n) + l1*(-3.0 + 4.0*n));
      ACoeff_t den = (-1.0*n + 2.0*n)*(l1 + 4.0*l1*n + 4.0*n*n);

      return this->normalizeDiag(n,-1)*(num/den);
   }

   D1Diags::ACoeff_t D1Diags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Ones(n.size());

      return this->normalizeDiag(n,0)*val;
   }

} // Chebyshev
} // Stencil
} // Worland
} // SparseSM
} // QuICC
