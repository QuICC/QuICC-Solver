/**
 * @file OpsQuadrature.hpp
 * @brief Mapping Generic Worland operator builders to specific Ops
 */

#pragma once

// External includes
//

// Project includes
//
#include "QuICC/Polynomial/Quadrature/WorlandRule.hpp"
#include "ViewOps/Worland/Tags.hpp"

namespace QuICC {
namespace Transform {
namespace Worland {

/// @brief This namespace hides implementation details
namespace details {

/// @brief Default Jones-Worland quadrature
/// @tparam TAG kind
template <class TAG>
struct OpsQuadratureMap {
    using type = ::QuICC::Polynomial::Quadrature::WorlandRule;
};

} // namespace details

/// @brief Convenience wrapper
/// (you cannot specialize type aliases)
/// @tparam TAG kind
template <class TAG>
using OpsQuadrature = typename details::OpsQuadratureMap<TAG>::type;

// // Specialized mapping for reductors
// namespace details {

// using namespace QuICC::Polynomial::Worland;


// } // namespace details

} // namespace Worland
} // namespace Transform
} // namespace QuICC
