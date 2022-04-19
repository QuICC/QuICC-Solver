/** 
 * @file I2Diags.cpp
 * @brief Source of the implementation of the full sphere Worland I2 sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// External includes
//

// Class include
//
#include "QuICC/SparseSM/Worland/I2Diags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I2Diags::I2Diags(const Scalar_t alpha, const Scalar_t dBeta, const int l)
      : IDiags(alpha, dBeta, l)
   {
   }

   I2Diags::~I2Diags()
   {
   }

}
}
}
