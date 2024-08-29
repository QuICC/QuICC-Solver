/**
 * @file Tags.hpp
 * @brief Tag types
 */
#pragma once

// System includes
//

// Project includes
//

namespace QuICC {
namespace Transform {
namespace ALegendre {

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

/// @brief Ll op type tag
struct Ll_t
{
};

/// @brief LlD1 op type tag
struct LlD1_t
{
};

/// @brief DivS1 op type tag
struct DivS1_t
{
};

/// @brief DivS1Dp op type tag
struct DivS1Dp_t
{
};

/// @brief LlDivS1 op type tag
struct LlDivS1_t
{
};

/// @brief LlDivS1Dp op type tag
struct LlDivS1Dp_t
{
};

/// @brief Ll2 op type tag
struct Ll2_t
{
};

/// @brief DivLl op type tag
struct DivLl_t
{
};

/// @brief DivLlD1 op type tag
struct DivLlD1_t
{
};

/// @brief DivLlDivS1 op type tag
struct DivLlDivS1_t
{
};

/// @brief DivLlDivS1Dp op type tag
struct DivLlDivS1Dp_t
{
};

} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
