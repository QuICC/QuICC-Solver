/** 
 * @file InsulatingSphereDiags.cpp
 * @brief Source of the implementation of the full sphere Worland insulating sphere boundary condition stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/Stencil/InsulatingSphereDiags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   InsulatingSphereDiags::InsulatingSphereDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l, 0)
   {
   }

} // Stencil
} // Worland
} // SparseSM
} // QuICC
