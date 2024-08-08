/**
 * @file Tags.hpp
 * @brief Tag types
 */
#pragma once

// External includes
//
#include <cstdint>

// Project includes
//

namespace QuICC {
namespace Transform {
namespace Worland {

//
// Tags
//

/// @brief tag type for projection direction.
/// Forwards i.e. physical to modal (integrator)
struct fwd_t
{
};

/// @brief tag type for projection direction.
/// Backwards i.e. modal to physical (projector)
struct bwd_t
{
};

/// @brief P op type tag
struct P_t
{
};

/// @brief D1 op type tag
struct D1_t
{
};

/// @brief D1R1 op type tag
struct D1R1_t
{
};


} // namespace Worland
} // namespace Transform
} // namespace QuICC
