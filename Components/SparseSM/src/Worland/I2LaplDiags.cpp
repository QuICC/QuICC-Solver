/** 
 * @file I2LaplDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I2Lapl sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Class include
//
#include "QuICC/SparseSM/Worland/I2LaplDiags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I2LaplDiags::I2LaplDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
