/**
 * @file Tags.hpp
 * @brief Tag types
 */
#pragma once

// System includes
//
#include <cstdint>
#include <array>

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
   constexpr static std::array<std::size_t,1> p = {0};
   constexpr static std::array<std::size_t,1> t = {0};
};

/// @brief tag type for spectral padding operation
struct spec_y1
{
   constexpr static std::array<std::size_t,1> p = {0};
   constexpr static std::array<std::size_t,1> t = {1};
};

/// @brief tag type for spectral padding operation
struct spec_d1
{
   constexpr static std::array<std::size_t,1> p = {1};
   constexpr static std::array<std::size_t,1> t = {0};
};

/// @brief tag type for spectral padding operation
struct spec_d2
{
   constexpr static std::array<std::size_t,1> p = {2};
   constexpr static std::array<std::size_t,1> t = {0};
};

/// @brief tag type for spectral padding operation
struct spec_d3
{
   constexpr static std::array<std::size_t,1> p = {3};
   constexpr static std::array<std::size_t,1> t = {0};
};

/// @brief tag type for spectral padding operation
struct spec_d4
{
   constexpr static std::array<std::size_t,1> p = {4};
   constexpr static std::array<std::size_t,1> t = {0};
};

/// @brief tag type for spectral padding operation
struct spec_d1y1
{
   constexpr static std::array<std::size_t,1> p = {1};
   constexpr static std::array<std::size_t,1> t = {1};
};

/// @brief tag type for spectral padding operation
struct spec_d1y2d1
{
   constexpr static std::array<std::size_t,2> p = {1,1};
   constexpr static std::array<std::size_t,2> t = {0,2};
};

/// @brief tag type for identity grid operation
struct grid_id
{
};

/// @brief tag type for division by Y^1 grid operation
struct grid_divy1
{
};

/// @brief tag type for division by Y^1 grid operation
struct grid_divy2
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
constexpr std::uint16_t none_t = 0;

/// @brief no special treatment
constexpr std::uint16_t ndealias_in = 1 << 1;

/// @brief no special treatment
constexpr std::uint16_t ndealias_out = 1 << 2;

/// @brief no special treatment
constexpr std::uint16_t zero_pad = 1 << 3;

/// @brief P op type tag
struct P_t
{
};

/// @brief P_Zero op type tag
struct P_Zero_t
{
};

/// @brief P op type tag
struct Y1_t
{
};

/// @brief DivY1_t op type tag
struct DivY1_t
{
};

/// @brief DivY2_t op type tag
struct DivY2_t
{
};

/// @brief D1_t op type tag
struct D1_t
{
};

/// @brief D2_t op type tag
struct D2_t
{
};

/// @brief D3_t op type tag
struct D3_t
{
};

/// @brief D4_t op type tag
struct D4_t
{
};

/// @brief D1Y1_t op type tag
struct D1Y1_t
{
};

/// @brief DivY1D1Y1_t op type tag
struct DivY1D1Y1_t
{
};

/// @brief SphRadLapl_t op type tag
struct SphRadLapl_t
{
};

} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
