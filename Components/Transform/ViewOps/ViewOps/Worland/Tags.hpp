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

/// @brief mask for special treatment.
/// none
constexpr std::uint16_t none_m = 0;

/// @brief mask for special treatment.
/// diff phi
constexpr std::uint16_t diffPhi_m = 1;


} // namespace Worland
} // namespace Transform
} // namespace QuICC
