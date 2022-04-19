/** 
 * @file R2.cpp
 * @brief Source of the implementation of the full sphere Worland R2 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/R2Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

   R2Diags::R2Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l)
   {
   }

   R2Diags::~R2Diags()
   {
   }

}
}
}
