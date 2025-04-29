/**
 * @file Types.hpp
 * @brief Fourier Mixed types
 */
#pragma once

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
namespace Transform {
namespace Fourier {
namespace Mixed {

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View::View<std::complex<double>, View::DCCSC3D>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View::View<double, View::DCCSC3D>;


} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
