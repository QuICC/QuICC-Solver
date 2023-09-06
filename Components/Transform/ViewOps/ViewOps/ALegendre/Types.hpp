/**
 * @file Types.hpp
 * @brief Associate Legendre types
 */
#pragma once

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
namespace Transform {
namespace ALegendre {

using namespace Memory;

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View<std::complex<double>, TRCLCSC3D>;
/// @brief mode coefficients view with row major layout (in a layer)
using modsRM_t = View<std::complex<double>, TRCLCSC3DJIK>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View<std::complex<double>, DCCSC3D>;
/// @brief physical coefficients view with row major layout (in a layer)
using physRM_t = View<std::complex<double>, DCCSC3DJIK>;
/// @brief projector view with column major layout (in a layer)
using proj_t = View<double, CTRRL3D>;
/// @brief projector view with row major layout (in a layer)
using projRM_t = View<double, CTRRL3DJIK>;
/// @brief projector view with column major layout (in a layer)
using int_t = View<double, TRCLCSC3D>;
/// @brief projector view with row major layout (in a layer)
using intRM_t = View<double, TRCLCSC3DJIK>;

} // namespace ALegendre
} // namespace Transform
} // namespace QuICC