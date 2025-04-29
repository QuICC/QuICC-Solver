/** 
 * @file InsulatingSphereDiags.cpp
 * @brief Source of the implementation of the full sphere Worland insulating sphereDiags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/InsulatingSphereDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace Chebyshev {

   InsulatingSphereDiags::InsulatingSphereDiags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::Stencil::InsulatingSphereDiags(alpha, MHD_MP(-0.5), l)
   {
   }

   InsulatingSphereDiags::ACoeff_t InsulatingSphereDiags::d_1(const ACoeff_t& n) const
   {
      auto l1 = this->l();

      ACoeff_t num = -2.0*n*(4.0*(-1.0 + n)*(-1.0 + n) + l1*(-2.0 + 4.0*n) + 1.0);
      ACoeff_t den = (-1.0 + 2.0*n)*(2.0*l1 + 1.0 + 4.0*l1*n + 4.0*n*n);

      ACoeff_t val = num/den;

      return this->normalizeDiag(n,-1)*val;
   }

   InsulatingSphereDiags::ACoeff_t InsulatingSphereDiags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Ones(n.size());

      return this->normalizeDiag(n,0)*val;
   }

} // Chebyshev
} // Stencil
} // Worland
} // SparseSM
} // QuICC
