/** 
 * @file ValueDiags.cpp
 * @brief Source of the implementation of the full sphere Worland ValueDiags sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Stencil/Chebyshev/ValueDiags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

namespace Chebyshev {

   ValueDiags::ValueDiags(const Scalar_t alpha, const int l)
      : QuICC::SparseSM::Worland::Stencil::ValueDiags(alpha, MHD_MP(-0.5), l)
   {
   }

   ValueDiags::ACoeff_t ValueDiags::d_1(const ACoeff_t& n) const
   {
      ACoeff_t val;

      val = -2.0*n/(2.0*n - 1.0);

      return this->normalizeDiag(n,-1)*val;
   }

   ValueDiags::ACoeff_t ValueDiags::d0(const ACoeff_t& n) const
   {
      ACoeff_t val = ACoeff_t::Ones(n.size());

      return this->normalizeDiag(n,0)*val;
   }

} // Chebyshev
} // Stencil
} // Worland
} // SparseSM
} // QuICC
