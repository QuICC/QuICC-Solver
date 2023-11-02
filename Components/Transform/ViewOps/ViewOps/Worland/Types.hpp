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
namespace Worland {

/// @brief namespace for uniform truncation
namespace Uniform {

using namespace Memory;

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View<std::complex<double>, DCCSC3D>;
/// @brief mode coefficients view with row major layout (in a layer)
using modsRM_t = View<std::complex<double>, DCCSC3DJIK>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View<std::complex<double>, DCCSC3D>;
/// @brief physical coefficients view with row major layout (in a layer)
using physRM_t = View<std::complex<double>, DCCSC3DJIK>;
/// @brief projector view with column major layout (in a layer)
using proj_t = View<double, CSL3D>;
/// @brief projector view with row major layout (in a layer)
using projRM_t = View<double, CSL3DJIK>;
/// @brief projector view with column major layout (in a layer)
using int_t = View<double, CSL3D>;
/// @brief projector view with row major layout (in a layer)
using intRM_t = View<double, CSL3DJIK>;

} // namespace Uniform

/// @brief namespace for triangular/trapezoidal truncation
namespace Triangular {

using namespace Memory;

/// @brief mode coefficients view with column major layout (in a layer)
using mods_t = View<std::complex<double>, S2CLCSC3D>;
/// @brief mode coefficients view with row major layout (in a layer)
using modsRM_t = View<std::complex<double>, S2CLCSC3DJIK>;
/// @brief physical coefficients view with column major layout (in a layer)
using phys_t = View<std::complex<double>, DCCSC3D>;
/// @brief physical coefficients view with row major layout (in a layer)
using physRM_t = View<std::complex<double>, DCCSC3DJIK>;
/// @brief projector view with column major layout (in a layer)
using proj_t = View<double, CS2RL3D>;
/// @brief projector view with row major layout (in a layer)
using projRM_t = View<double, CS2RL3DJIK>;
/// @brief projector view with column major layout (in a layer)
using int_t = View<double, S2CLCSC3D>;
/// @brief projector view with row major layout (in a layer)
using intRM_t = View<double, S2CLCSC3DJIK>;

} // namespace Triangular

} // namespace Worland
} // namespace Transform
} // namespace QuICC
