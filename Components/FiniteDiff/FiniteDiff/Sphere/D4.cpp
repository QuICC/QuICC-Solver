/**
 * @file D4.cpp
 * @brief Source of D4 operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/D4.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   D4::D4(const std::size_t order, const std::size_t zTop, const std::size_t zBot)
      : Operator(order, zTop, zBot)
   {
   }

   D4::D4(const std::size_t zTop, const std::size_t zBot)
      : Operator(QUICC_FINITEDIFF_SPHERE_ORDER, zTop, zBot)
   {
   }

   D4::D4()
      : D4(QUICC_FINITEDIFF_SPHERE_ORDER, 0, 0)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
