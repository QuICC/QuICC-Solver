/**
 * @file Tags.hpp
 * @brief Tag types
 */
#pragma once

// System includes
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

} // namespace Worland
} // namespace Transform
} // namespace QuICC
