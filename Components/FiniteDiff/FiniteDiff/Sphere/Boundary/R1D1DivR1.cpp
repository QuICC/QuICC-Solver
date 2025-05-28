/**
 * @file R1D1DivR1.cpp
 * @brief Source of R1D1DivR1 boundary operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Boundary/R1D1DivR1.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

namespace Boundary {

   R1D1DivR1::R1D1DivR1(const std::size_t order)
      : Operator(order)
   {
   }

   R1D1DivR1::R1D1DivR1()
      : R1D1DivR1(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Boundary
} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
