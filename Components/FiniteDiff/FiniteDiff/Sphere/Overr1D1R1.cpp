/**
 * @file Overr1D1R1.cpp
 * @brief Source of Overr1D1R1 operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Overr1D1R1.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   Overr1D1R1::Overr1D1R1(const std::size_t order)
      : Operator(order)
   {
   }

   Overr1D1R1::Overr1D1R1()
      : Overr1D1R1(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
