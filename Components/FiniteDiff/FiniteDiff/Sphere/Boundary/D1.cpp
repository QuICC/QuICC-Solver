/**
 * @file D1.cpp
 * @brief Source of D1 boundary operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Boundary/D1.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

namespace Boundary {

   D1::D1(const std::size_t order)
      : Operator(order)
   {
   }

   D1::D1()
      : D1(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Boundary
} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
