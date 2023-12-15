/**
 * @file Types.hpp
 * @brief Fourier Complex types
 */
#pragma once

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Complex {

/// @brief mode coefficients view with column major layout (in a layer)
/// with in-order mapping
using mods_t = View::View<std::complex<double>, View::DCCSC3DInOrder>;
/// @brief physical coefficients view with column major layout (in a layer)
/// with in-order mapping
using phys_t = View::View<std::complex<double>, View::DCCSC3DInOrder>;

} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
