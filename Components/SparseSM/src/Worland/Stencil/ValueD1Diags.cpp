/** 
 * @file ValueD1Diags.cpp
 * @brief Source of the implementation of the full sphere Worland ValueD1 boundary condition stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/ValueD1Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   ValueD1Diags::ValueD1Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l, 0)
   {
   }

} // Stencil
} // Worland
} // SparseSM
} // QuICC
