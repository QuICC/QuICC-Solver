/**
 * @file Types.hpp
 * @brief Worland types
 */
#pragma once

// Project includes
//
#include "View/View.hpp"

namespace QuICC {
namespace Transform {
/// @brief namespace for Worland operators
namespace Worland {
/// @brief namespace for uniform truncation
namespace Uniform {

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View::View<std::complex<double>, View::DCCSC3D>;
/// @brief mode coefficients view with row major layout (in a layer)
using modsRM_t = View::View<std::complex<double>, View::DCCSC3DJIK>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View::View<std::complex<double>, View::DCCSC3D>;
/// @brief physical coefficients view with row major layout (in a layer)
using physRM_t = View::View<std::complex<double>, View::DCCSC3DJIK>;
/// @brief projector view with column major layout (in a layer)
using proj_t = View::View<double, View::CSL3D>;
/// @brief projector view with row major layout (in a layer)
using projRM_t = View::View<double, View::CSL3DJIK>;
/// @brief projector view with column major layout (in a layer)
using int_t = View::View<double, View::CSL3D>;
/// @brief projector view with row major layout (in a layer)
using intRM_t = View::View<double, View::CSL3DJIK>;

} // namespace Uniform

/// @brief namespace for triangular/trapezoidal truncation
namespace Triangular {

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View::View<std::complex<double>, View::S2CLCSC3D>;
/// @brief mode coefficients view with row major layout (in a layer)
using modsRM_t = View::View<std::complex<double>, View::S2CLCSC3DJIK>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View::View<std::complex<double>, View::DCCSC3D>;
/// @brief physical coefficients view with row major layout (in a layer)
using physRM_t = View::View<std::complex<double>, View::DCCSC3DJIK>;
/// @brief projector view with column major layout (in a layer)
using proj_t = View::View<double, View::CS2RL3D>;
/// @brief projector view with row major layout (in a layer)
using projRM_t = View::View<double, View::CS2RL3DJIK>;
/// @brief projector view with column major layout (in a layer)
using int_t = View::View<double, View::S2CLCSC3D>;
/// @brief projector view with row major layout (in a layer)
using intRM_t = View::View<double, View::S2CLCSC3DJIK>;

} // namespace Triangular

} // namespace Worland
} // namespace Transform
} // namespace QuICC
