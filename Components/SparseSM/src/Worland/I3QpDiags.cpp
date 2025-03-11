/** 
 * @file I3QpDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I3Qp sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I3QpDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I3QpDiags::I3QpDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
