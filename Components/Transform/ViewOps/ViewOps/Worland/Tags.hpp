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

/// @brief P_Zero op type tag
struct P_Zero_t
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

/// @brief R1_Zero op type tag
struct R1_Zero_t
{
};

/// @brief I2 op type tag
struct I2_t
{
};

/// @brief I2_Zero op type tag
struct I2_Zero_t
{
};

/// @brief I2DivR1_Zero op type tag
struct I2DivR1_Zero_t
{
};

/// @brief I4DivR1_Zero op type tag
struct I4DivR1_Zero_t
{
};

/// @brief I6DivR1_Zero op type tag
struct I6DivR1_Zero_t
{
};

/// @brief I2DivR1D1R1_Zero op type tag
struct I2DivR1D1R1_Zero_t
{
};

/// @brief I4DivR1D1R1_Zero op type tag
struct I4DivR1D1R1_Zero_t
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


} // namespace Worland
} // namespace Transform
} // namespace QuICC
