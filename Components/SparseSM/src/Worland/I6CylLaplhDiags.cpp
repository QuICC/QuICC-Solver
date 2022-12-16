/** 
 * @file I6CylLaplhDiags.cpp
 * @brief Source of the implementation of the full sphere Worland I6CylLaplh sparse operator
 */

// System includes
//
#include <cassert>
#include <stdexcept>

// Class include
//
#include "QuICC/SparseSM/Worland/I6CylLaplhDiags.hpp"

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Worland {

   I6CylLaplhDiags::I6CylLaplhDiags(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, q)
   {
   }

}
}
}
