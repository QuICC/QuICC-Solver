/**
 * @file OpsBuilder.hpp
 * @brief Mapping Generic Worland operator builders to specific Ops
 */

#pragma once

// External includes
//

// Project includes
//
#include "QuICC/Polynomial/Worland/Wnl.hpp"
#include "QuICC/Polynomial/Worland/dWnl.hpp"
#include "QuICC/Polynomial/Worland/drWnl.hpp"
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Worland/Builder.hpp"

namespace QuICC {
namespace Transform {
namespace Worland {

/// @brief This namespace hides implementation details
namespace details {

/// @brief Generic mapping
/// @tparam VOP operator view type
/// @tparam TAG kind
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class TAG, class DIR>
struct OpsBuilderMap {
    using type = void;
};

} // namespace details

/// @brief Convenience wrapper
/// (you cannot specialize type aliases)
/// @tparam VOP operator view type
/// @tparam TAG kind
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class TAG, class DIR>
using OpsBuilder = typename details::OpsBuilderMap<VOP, TAG, DIR>::type;

// Actual mapping
namespace details {

using namespace QuICC::Polynomial::Worland;

/// @brief P Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, P_t, DIR> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<Wnl>, DIR>;
};

/// @brief D1 Builder
/// @tparam VOP operator view type
/// @tparam DIR fwd_t or bwd_t
template <class VOP, class DIR>
struct OpsBuilderMap<VOP, D1_t, DIR> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<dWnl>, DIR>;
};

/// @brief D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, D1R1_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<drWnl>, bwd_t>;
};

} // namespace details


} // namespace Worland
} // namespace Transform
} // namespace QuICC
