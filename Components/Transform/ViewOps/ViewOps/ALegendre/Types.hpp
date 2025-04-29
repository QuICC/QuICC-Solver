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
/// @brief namespace for Associated Legendre operators
namespace ALegendre {

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View::View<std::complex<double>, View::S1CLCSC3D>;
/// @brief mode coefficients view with row major layout (in a layer)
using modsRM_t = View::View<std::complex<double>, View::S1CLCSC3DJIK>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View::View<std::complex<double>, View::DCCSC3D>;
/// @brief physical coefficients view with row major layout (in a layer)
using physRM_t = View::View<std::complex<double>, View::DCCSC3DJIK>;
/// @brief projector view with column major layout (in a layer)
using proj_t = View::View<double, View::CS1RL3D>;
/// @brief projector view with row major layout (in a layer)
using projRM_t = View::View<double, View::CS1RL3DJIK>;
/// @brief integrator view with column major layout (in a layer)
using int_t = View::View<double, View::S1CLCSC3D>;
/// @brief integrator view with row major layout (in a layer)
using intRM_t = View::View<double, View::S1CLCSC3DJIK>;

} // namespace ALegendre
} // namespace Transform
} // namespace QuICC
