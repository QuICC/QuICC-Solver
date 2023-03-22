/** 
 * @file ValueDiags.cpp
 * @brief Source of the implementation of the full sphere Worland Value boundary condition stencil
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/Stencil/ValueDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Stencil {

   ValueDiags::ValueDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l, 0)
   {
   }

} // Stencil
} // Worland
} // SparseSM
} // QuICC
