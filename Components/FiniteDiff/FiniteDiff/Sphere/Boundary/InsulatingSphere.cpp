/**
 * @file InsulatingSphere.cpp
 * @brief Source of InsulatingSphere boundary operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Boundary/InsulatingSphere.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

namespace Boundary {

   InsulatingSphere::InsulatingSphere(const std::size_t order)
      : Operator(order)
   {
   }

   InsulatingSphere::InsulatingSphere()
      : InsulatingSphere(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Boundary
} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
