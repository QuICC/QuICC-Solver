/** 
 * @file I4Lapl2Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I4Lapl2 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4Lapl2Diags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I4Lapl2Diags::I4Lapl2Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
