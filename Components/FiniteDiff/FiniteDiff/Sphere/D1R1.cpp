/**
 * @file D1R1.cpp
 * @brief Source of D1R1 operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/D1R1.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   D1R1::D1R1(const std::size_t order, const std::size_t zTop, const std::size_t zBot)
      : Operator(order, zTop, zBot)
   {
   }

   D1R1::D1R1(const std::size_t zTop, const std::size_t zBot)
      : D1R1(QUICC_FINITEDIFF_SPHERE_ORDER, zTop, zBot)
   {
   }

   D1R1::D1R1()
      : D1R1(QUICC_FINITEDIFF_SPHERE_ORDER, 1, 1)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
