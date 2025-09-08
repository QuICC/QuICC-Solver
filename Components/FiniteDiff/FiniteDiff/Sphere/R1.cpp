/**
 * @file R1.cpp
 * @brief Source of R1 operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/R1.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   R1::R1(const std::size_t order, const std::size_t zTop, const std::size_t zBot)
      : Operator(order, zTop, zBot)
   {
   }

   R1::R1(const std::size_t zTop, const std::size_t zBot)
      : Operator(QUICC_FINITEDIFF_SPHERE_ORDER, zTop, zBot)
   {
   }

   R1::R1()
      : R1(QUICC_FINITEDIFF_SPHERE_ORDER, 0, 0)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
