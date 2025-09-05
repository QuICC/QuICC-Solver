/**
 * @file D3.cpp
 * @brief Source of D3 operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/D3.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   D3::D3(const std::size_t order)
      : Operator(order)
   {
   }

   D3::D3()
      : D3(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
