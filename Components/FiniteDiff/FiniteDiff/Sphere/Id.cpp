/**
 * @file Id.cpp
 * @brief Source of Id operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/Id.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   Id::Id(const std::size_t order, const std::size_t zTop, const std::size_t zBot)
      : Operator(order, zTop, zBot)
   {
   }

   Id::Id(const std::size_t zTop, const std::size_t zBot)
      : Id(QUICC_FINITEDIFF_SPHERE_ORDER, zTop, zBot)
   {
   }

   Id::Id()
      : Id(QUICC_FINITEDIFF_SPHERE_ORDER, 0, 0)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
