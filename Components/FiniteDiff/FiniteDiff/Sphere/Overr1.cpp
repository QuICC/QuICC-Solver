/**
 * @file Overr1.cpp
 * @brief Source of Overr1 operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Overr1.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   Overr1::Overr1(const std::size_t order)
      : Operator(order)
   {
   }

   Overr1::Overr1()
      : Overr1(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
