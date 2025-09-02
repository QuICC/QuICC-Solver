/**
 * @file IntegralR2.cpp
 * @brief Source of Integral r^2 operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/IntegralR2.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   IntegralR2::IntegralR2(const std::size_t order)
      : Operator(order), mcTaylorN(4)
   {
   }

   IntegralR2::IntegralR2()
      : IntegralR2(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
