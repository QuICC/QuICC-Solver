/** 
 * @file ValueD2Diags.cpp
 * @brief Source of the implementation of the full sphere Worland ValueD2 boundary condition stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/ValueD2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   ValueD2Diags::ValueD2Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l, 0)
   {
   }

} // Stencil
} // Worland
} // SparseSM
} // QuICC
