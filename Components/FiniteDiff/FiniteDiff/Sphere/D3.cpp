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

   D3::D3(const std::size_t order, const std::size_t zTop, const std::size_t zBot)
      : Operator(order, zTop, zBot)
   {
   }

   D3::D3(const std::size_t zTop, const std::size_t zBot)
      : Operator(QUICC_FINITEDIFF_SPHERE_ORDER, zTop, zBot)
   {
   }

   D3::D3()
      : D3(QUICC_FINITEDIFF_SPHERE_ORDER, 0, 0)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
