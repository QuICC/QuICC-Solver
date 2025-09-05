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

   Overr1::Overr1(const std::size_t order, const std::size_t zTop, const std::size_t zBot)
      : Operator(order, zTop, zBot)
   {
   }

   Overr1::Overr1(const std::size_t zTop, const std::size_t zBot)
      : Operator(QUICC_FINITEDIFF_SPHERE_ORDER, zTop, zBot)
   {
   }

   Overr1::Overr1()
      : Overr1(QUICC_FINITEDIFF_SPHERE_ORDER, 1, 1)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
