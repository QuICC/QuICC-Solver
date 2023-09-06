/** 
 * @file I4Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I4 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I4Diags::I4Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
