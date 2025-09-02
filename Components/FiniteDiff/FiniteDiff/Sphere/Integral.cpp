/**
 * @file Integral.cpp
 * @brief Source of Integral operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Integral.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   Integral::Integral(const std::size_t order)
      : Operator(order), mcTaylorN(4)
   {
   }

   Integral::Integral()
      : Integral(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
