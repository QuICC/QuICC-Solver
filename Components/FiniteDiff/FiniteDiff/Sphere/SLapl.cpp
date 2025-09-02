/**
 * @file SLapl.cpp
 * @brief Source of SLapl operator
 */

// System includes
//

// Project includes
//
#include "FiniteDiff/Sphere/SLapl.hpp"

namespace QuICC {

namespace FiniteDiff {

namespace Sphere {

   SLapl::SLapl(const std::size_t order, const std::size_t zTop, const std::size_t zBot)
      : Operator(order, zTop, zBot)
   {
   }

   SLapl::SLapl(const std::size_t zTop, const std::size_t zBot)
      : SLapl(QUICC_FINITEDIFF_SPHERE_ORDER, zTop, zBot)
   {
   }

   SLapl::SLapl()
      : SLapl(QUICC_FINITEDIFF_SPHERE_ORDER, 1, 1)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
