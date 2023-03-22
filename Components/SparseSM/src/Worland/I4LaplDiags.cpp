/** 
 * @file I4LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I4Lapl sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I4LaplDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I4LaplDiags::I4LaplDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
