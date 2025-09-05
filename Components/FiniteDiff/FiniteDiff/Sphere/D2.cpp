/**
 * @file D2.cpp
 * @brief Source of D2 operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/D2.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   D2::D2(const std::size_t order)
      : Operator(order)
   {
   }

   D2::D2()
      : D2(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
