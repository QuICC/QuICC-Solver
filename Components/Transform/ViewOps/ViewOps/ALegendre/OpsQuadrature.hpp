/**
 * @file OpsQuadrature.hpp
 * @brief Mapping Generic ALegendre operator builders to specific Ops
 */

#pragma once

// External includes
//

// Project includes
//
#include "QuICC/Polynomial/Quadrature/LegendreRule.hpp"
#include "Types/Internal/Math.hpp"
#include "ViewOps/ALegendre/Tags.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {


/// @brief This namespace hides implementation details
namespace details {

struct ALegendreRule
{
   void computeQuadrature(Internal::Array& igrid, Internal::Array& iweights,
      const int gSize);
};


/// @brief Default ALegendre quadrature
/// @tparam TAG kind
template <class TAG> struct OpsQuadratureMap
{
   using type = ALegendreRule;
};

} // namespace details

/// @brief Convenience wrapper
/// (you cannot specialize type aliases)
/// @tparam TAG kind
template <class TAG>
using OpsQuadrature = typename details::OpsQuadratureMap<TAG>::type;


} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
