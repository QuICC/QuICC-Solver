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
namespace Chebyshev {
namespace LinearMap {

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

/// @brief tag type for identity spectral operation
struct spec_id
{
};

/// @brief tag type for identity grid operation
struct grid_id
{
};

///view cpu implementation tag
struct viewCpu_t
{
};

/// view gpu implementation tag
struct viewGpu_t
{
};

/// view gpu VkFFT implementation tag
struct viewGpuVkFFT_t
{
};

/// @brief no special treatment
constexpr std::uint16_t none_l = 0;

/// @brief P op type tag
struct P_t
{
};

/// @brief P_Zero op type tag
struct P_Zero_t
{
};


} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
