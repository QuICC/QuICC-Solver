/** 
 * @file ValueD1Diags.cpp
 * @brief Source of the implementation of the full sphere Worland ValueD1 boundary condition stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Stencil/ValueD1Diags.hpp"

// Project includes
//

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
