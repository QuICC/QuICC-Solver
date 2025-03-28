/**
 * @file Types.hpp
 * @brief Chebyshev LinearMap types
 */
#pragma once

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
namespace Transform {
namespace Chebyshev {
namespace LinearMap {

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View::View<std::complex<double>, View::DCCSC3D>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View::View<std::complex<double>, View::DCCSC3D>;


} // namespace LinearMap
} // namespace Chebyshev
} // namespace Transform
} // namespace QuICC
