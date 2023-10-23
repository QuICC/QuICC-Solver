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

using namespace Memory;

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View<std::complex<double>, DCCSC3D>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View<double, DCCSC3D>;


} // namespace Mixed
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
