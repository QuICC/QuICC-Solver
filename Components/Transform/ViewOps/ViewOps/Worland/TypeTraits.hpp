/**
 * @file TypeTraits.hpp
 * @brief Worland types
 */
#pragma once

// Project includes
//
#include "ViewOps/Worland/Tags.hpp"
#include "ViewOps/Worland/Types.hpp"

namespace QuICC {
namespace Transform {
namespace Worland {

namespace Uniform {

/// @brief Generic template to check if an operator is a projector
/// @tparam TView
/// @tparam TDirection
template <class TView, class TDirection> struct is_projector : std::false_type
{
};

/// @brief Specialization for row major projector
template <> struct is_projector<projRM_t, bwd_t> : std::true_type
{
};

/// @brief Specialization for col major projector
template <> struct is_projector<proj_t, bwd_t> : std::true_type
{
};

/// @brief Helper
/// @tparam TView
/// @tparam TDirection
template <class TView, class TDirection>
inline constexpr bool is_projector_v = is_projector<TView, TDirection>::value;

/// @brief Generic template to check if an operator is an integrator
/// @tparam TView
/// @tparam TDirection
template <class TView, class TDirection> struct is_integrator : std::false_type
{
};

/// @brief Specialization for row major integrator
template <> struct is_integrator<intRM_t, fwd_t> : std::true_type
{
};

/// @brief Specialization for col major integrator
template <> struct is_integrator<int_t, fwd_t> : std::true_type
{
};

/// @brief Helper
/// @tparam TView
/// @tparam TDirection
template <class TView, class TDirection>
inline constexpr bool is_integrator_v = is_integrator<TView, TDirection>::value;


} // namespace Uniform

} // namespace Worland
} // namespace Transform
} // namespace QuICC
