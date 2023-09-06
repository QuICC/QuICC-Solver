/** 
 * @file I4QmDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4Qm sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4QmDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I4QmDiags::I4QmDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
