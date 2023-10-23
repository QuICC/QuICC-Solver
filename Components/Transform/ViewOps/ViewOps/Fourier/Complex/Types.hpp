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

using namespace Memory;

/// @brief mode coefficients view with column major layout (in a layer)
/// with in-order mapping
using mods_t = View<std::complex<double>, DCCSC3DInOrder>;
/// @brief physical coefficients view with column major layout (in a layer)
/// with in-order mapping
using phys_t = View<std::complex<double>, DCCSC3DInOrder>;

} // namespace Complex
} // namespace Fourier
} // namespace Transform
} // namespace QuICC
