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
namespace Fourier {

//
// Tags
//

/// @brief tag type for differentiation direction.
/// Forwards i.e. physical to modal
struct fwd_t
{
};

/// @brief tag type for differentiation direction.
/// Backwards i.e. modal to physical
struct bwd_t
{
};

/// @brief mask for special treatment.
/// none
constexpr std::uint16_t none_m = 0;

/// @brief mask for special treatment.
/// zero mode is always zero just projection
constexpr std::uint16_t zeroP_m = 1 << 1;

/// @brief mask for special treatment.
/// zero mode is always zero just projection
/// with negated sign
constexpr std::uint16_t zeroMinusP_m = 1 << 2;

/// @brief mask for special treatment.
/// m = 0, k = 0 n != 0 rest to zero,
/// i.e. keep mean only
constexpr std::uint16_t zeroResetMean_m = 1 << 3;

/// @brief mask for special treatment.
/// inverse (1/coeff)
constexpr std::uint16_t inverse_m = 1 << 4;


} // namespace Fourier
} // namespace Transform
} // namespace QuICC
