/** 
 * @file I4Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I4 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/I4Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I4Diags::I4Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l)
   {
   }

   I4Diags::~I4Diags()
   {
   }

}
}
}
