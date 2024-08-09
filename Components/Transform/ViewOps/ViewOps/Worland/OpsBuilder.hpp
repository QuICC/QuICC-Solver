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
#include "QuICC/Polynomial/Worland/r_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/dr_1drWnl.hpp"
#include "QuICC/Polynomial/Worland/r_1Wnl.hpp"
#include "QuICC/Polynomial/Worland/slaplWnl.hpp"
#include "QuICC/Polynomial/Worland/claplhWnl.hpp"
#include "QuICC/Polynomial/Worland/r_1claplhWnl.hpp"
#include "QuICC/Polynomial/Worland/dclaplhWnl.hpp"
#include "DenseSM/Worland/Operator.hpp"
#include "DenseSM/Worland/OperatorWithMean.hpp"
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

/// @brief D1_P Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, D1_P_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::OperatorWithMean<dWnl, Wnl>, bwd_t>;
};

/// @brief DivR1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivR1_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<r_1Wnl<recurrence_t>>, bwd_t>;
};

/// @brief DivR1_Zero Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivR1_Zero_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::OperatorWithMean<r_1Wnl<recurrence_t>, void>, bwd_t>;
};

/// @brief DivR1D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivR1D1R1_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<r_1drWnl<recurrence_t>>, bwd_t>;
};

/// @brief DivR1D1R1_Zero Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivR1D1R1_Zero_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::OperatorWithMean<r_1drWnl<recurrence_t>, void>, bwd_t>;
};

/// @brief SphLapl Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, SphLapl_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<slaplWnl>, bwd_t>;
};

/// @brief CylLaplh Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, CylLaplh_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<claplhWnl>, bwd_t>;
};

/// @brief CylLaplh_DivR1D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, CylLaplh_DivR1D1R1_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::OperatorWithMean<claplhWnl, r_1drWnl<recurrence_t>>, bwd_t>;
};

/// @brief D1CylLaplh Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, D1CylLaplh_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::Operator<dclaplhWnl>, bwd_t>;
};

/// @brief D1CylLaplh_D1DivR1D1R1 Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, D1CylLaplh_D1DivR1D1R1_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::OperatorWithMean<dclaplhWnl, dr_1drWnl>, bwd_t>;
};

/// @brief DivR1CylLaplh_Zero Builder
/// Projector only
/// @tparam VOP operator view type
template <class VOP>
struct OpsBuilderMap<VOP, DivR1CylLaplh_Zero_t, bwd_t> {
    using type = typename Worland::Builder<VOP, QuICC::DenseSM::Worland::OperatorWithMean<r_1claplhWnl, void>, bwd_t>;
};

} // namespace details

} // namespace Worland
} // namespace Transform
} // namespace QuICC
