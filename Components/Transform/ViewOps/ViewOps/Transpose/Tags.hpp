/**
 * @file Tags.hpp
 * @brief Tag types
 */
#pragma once

// External includes
//

// Project includes
//

namespace QuICC {
namespace Transpose {

//
// Tags
//

/// @brief tag type for transpose permutation order.
/// perm = [0, 1, 2]
struct p012_t {};

/// @brief tag type for transpose permutation order.
/// perm = [2, 0, 1]
struct p201_t {};

/// @brief tag type for transpose permutation order.
/// perm = [1, 2, 0]
struct p120_t {};



} // namespace Transpose
} // namespace QuICC
