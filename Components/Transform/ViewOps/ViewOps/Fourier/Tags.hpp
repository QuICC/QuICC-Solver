/**
 * @file Tags.hpp
 * @brief Tag types
 */
#pragma once

// External includes
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
struct fwd_t {};

/// @brief tag type for differentiation direction.
/// Backwards i.e. modal to physical
struct bwd_t {};

/// @brief tag type for special treatment.
/// zero mode is always zero just projection
struct zeroP_t {};

/// @brief tag type for special treatment.
/// zero mode is always zero just projection
/// with negated sign
struct zeroMinusP_t {};

/// @brief tag type for special treatment.
/// m = 0, k = 0 n != 0 rest to zero,
/// i.e. keep mean only
struct zeroResetMean_t {};

/// @brief tag type for special treatment.
/// inverse (1/coeff)
struct inverse_t {};


} // namespace Fourier
} // namespace Transform
} // namespace QuICC
