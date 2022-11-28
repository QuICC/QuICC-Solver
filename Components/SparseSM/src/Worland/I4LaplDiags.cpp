/** 
 * @file I4LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4Lapl sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/I4LaplDiags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I4LaplDiags::I4LaplDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l)
   {
   }

}
}
}
