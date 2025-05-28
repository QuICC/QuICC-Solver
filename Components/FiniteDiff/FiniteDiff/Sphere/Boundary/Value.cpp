/**
 * @file Value.cpp
 * @brief Source of Value boundary operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Boundary/Value.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

namespace Boundary {

   Value::Value(const std::size_t order)
      : Operator(order)
   {
   }

   Value::Value()
      : Value(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Boundary
} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
