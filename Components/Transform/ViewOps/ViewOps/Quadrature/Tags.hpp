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
namespace Quadrature {

//
// Tags
//

/// @brief mask for special treatment.
/// none
constexpr std::uint16_t none_m = 0;

/// @brief mask for special treatment.
/// differentiate in phi for a projector
constexpr std::uint16_t diffPhiPrj_m = 1;

/// @brief mask for special treatment.
/// differentiate in phi for a integrator
constexpr std::uint16_t diffPhiInt_m = 1 << 1;


} // namespace Quadrature
} // namespace Transform
} // namespace QuICC
