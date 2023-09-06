/** 
 * @file I2QpDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I2Qp sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Project includes
//
#include "QuICC/SparseSM/Worland/I2QpDiags.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I2QpDiags::I2QpDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
