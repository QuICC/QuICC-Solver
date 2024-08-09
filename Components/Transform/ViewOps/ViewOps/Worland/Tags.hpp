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

/// @brief D1_P op type tag
struct D1_P_t
{
};

/// @brief DivR1 op type tag
struct DivR1_t
{
};

/// @brief DivR1_Zero op type tag
struct DivR1_Zero_t
{
};

/// @brief DivR1D1R1 op type tag
struct DivR1D1R1_t
{
};

/// @brief DivR1D1R1_Zero op type tag
struct DivR1D1R1_Zero_t
{
};

/// @brief SphLapl op type tag
struct SphLapl_t
{
};

/// @brief CylLaplh op type tag
struct CylLaplh_t
{
};

/// @brief CylLaplh_DivR1D1R1 op type tag
struct CylLaplh_DivR1D1R1_t
{
};

/// @brief D1CylLaplh op type tag
struct D1CylLaplh_t
{
};

/// @brief D1CylLaplh_D1DivR1D1R1 op type tag
struct D1CylLaplh_D1DivR1D1R1_t
{
};

/// @brief DivR1CylLaplh_Zero op type tag
struct DivR1CylLaplh_Zero_t
{
};

} // namespace Worland
} // namespace Transform
} // namespace QuICC
