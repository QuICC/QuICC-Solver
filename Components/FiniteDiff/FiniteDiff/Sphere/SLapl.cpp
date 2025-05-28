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

   SLapl::SLapl(const std::size_t order)
      : Operator(order)
   {
   }

   SLapl::SLapl()
      : SLapl(QUICC_FINITEDIFF_SPHERE_ORDER)
   {
   }

} // namespace Sphere
} // namespace FiniteDiff
} // namespace QuICC
