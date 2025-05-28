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

   Id::Id(const std::size_t order)
      : Operator(order)
   {
   }

   Id::Id()
      : Id(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
