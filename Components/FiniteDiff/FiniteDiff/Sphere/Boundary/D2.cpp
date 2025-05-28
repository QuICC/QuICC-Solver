/**
 * @file D2.cpp
 * @brief Source of D2 boundary operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Boundary/D2.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

namespace Boundary {

   D2::D2(const std::size_t order)
      : Operator(order)
   {
   }

   D2::D2()
      : D2(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Boundary
} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
