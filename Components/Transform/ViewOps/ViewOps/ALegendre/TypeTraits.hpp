/**
 * @file TypeTraits.hpp
 * @brief Associate Legendre types
 */
#pragma once

// Project includes
//
#include "ViewOps/ALegendre/Types.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {

/// @brief Generic template to check if an operator is a projector
/// @tparam T
template <class T>
struct is_projector : std::false_type {};

/// @brief Specialization for row major projector
template <>
struct is_projector<projRM_t> : std::true_type {};

/// @brief Specialization for col major projector
template <>
struct is_projector<proj_t> : std::true_type {};

/// @brief Helper
/// @tparam T
template <class T>
inline constexpr bool is_projector_v = is_projector<T>::value;

/// @brief Generic template to check if an operator is an integrator
/// @tparam T
template <class T>
struct is_integrator : std::false_type {};

/// @brief Specialization for row major integrator
template <>
struct is_integrator<intRM_t> : std::true_type {};

/// @brief Specialization for col major integrator
template <>
struct is_integrator<int_t> : std::true_type {};

/// @brief Helper
/// @tparam T
template <class T>
inline constexpr bool is_integrator_v = is_integrator<T>::value;


} // namespace ALegendre
} // namespace Transform
} // namespace QuICC