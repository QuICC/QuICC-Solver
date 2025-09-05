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
namespace FiniteDiff {
namespace Sphere {

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

/// @brief P_Zero op type tag
struct P_Zero_t
{
};

/// @brief D1 op type tag
struct D1_t
{
};

/// @brief D1 op type tag
struct D2_t
{
};

/// @brief D1 op type tag
struct D3_t
{
};

/// @brief D1 op type tag
struct D4_t
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

/// @brief R1 op type tag
struct R1_t
{
};

/// @brief R1_Zero op type tag
struct R1_Zero_t
{
};

/// @brief Energy op type tag
struct Energy_t
{
};

/// @brief EnergyD1R1 op type tag
struct EnergyD1R1_t
{
};

/// @brief EnergyR2 op type tag
struct EnergyR2_t
{
};

/// @brief EnergySLaplR2 op type tag
struct EnergySLaplR2_t
{
};

/// @brief Power op type tag
struct Power_t
{
};

/// @brief PowerD1R1 op type tag
struct PowerD1R1_t
{
};

/// @brief PowerR2 op type tag
struct PowerR2_t
{
};

/// @brief PowerSLaplR2 op type tag
struct PowerSLaplR2_t
{
};

/// @brief RadialPower op type tag
struct RadialPower_t
{
};

/// @brief RadialPowerDivR1 op type tag
struct RadialPowerDivR1_t
{
};

/// @brief RadialPowerDivR1D1R1 op type tag
struct RadialPowerDivR1D1R1_t
{
};


} // namespace Sphere
} // namespace FiniteDiff
} // namespace Transform
} // namespace QuICC
