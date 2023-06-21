/** 
 * @file Condition.cpp
 * @brief Source of the interface to a generic Worland boundary condition
 */

// System includes
//

// Project include
//
#include "QuICC/SparseSM/Worland/Boundary/ICondition.hpp"

namespace QuICC {

namespace SparseSM {

namespace Worland {

namespace Boundary {

   ICondition::ICondition(const Scalar_t alpha, const Scalar_t dBeta, const int l, const int q)
      : IDiags(alpha, dBeta, l, 0)
   {
   }
 
} // Boundary
} // Worland
} // Polynomial
} // QuICC
