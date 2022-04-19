/** 
 * @file I6Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I6 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/I6Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I6Diags::I6Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l)
   {
   }

   I6Diags::~I6Diags()
   {
   }

}
}
}
