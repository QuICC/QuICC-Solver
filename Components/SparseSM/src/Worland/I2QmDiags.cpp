/** 
 * @file I2QmDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I2Qm sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I2QmDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I2QmDiags::I2QmDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
