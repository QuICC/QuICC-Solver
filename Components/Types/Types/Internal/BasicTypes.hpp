/**
 * @file BasicTypes.hpp
 * @brief Some scalar typedefs for potentially higher precision
 * computation.
 * The switch between MP or not is a compile time option.
 */

#ifndef QUICC_TYPES_INTERNAL_BASICTYPES_HPP
#define QUICC_TYPES_INTERNAL_BASICTYPES_HPP

// System includes
//

// Project includes
//
#ifdef QUICC_MULTPRECISION
#include "Types/MP/Typedefs.hpp"
#else
#include "Types/Typedefs.hpp"
#endif

namespace QuICC {
/// @brief namespace for Interal types and operations
/// i.e. either QuICC or MP types, the switch is a compile time option
namespace Internal {

#ifdef QUICC_MULTPRECISION
/// @brief Is internal multiple precision?
constexpr bool isMP = true;
/// @brief Typedef for internal floating point type (MP)
using MHDFloat = MP::MHDFloat;
/// @brief Typedef for internal complex type value (MP)
using MHDComplex = MP::MHDComplex;
/// @brief Typedef for internal extended floating point type (MP)
using MHDLong = MP::MHDLong;
#else
/// @brief Is internal multiple precision?
constexpr bool isMP = false;
/// @brief Typedef for internal floating point type (non-MP)
using MHDFloat = MHDFloat;
/// @brief Typedef for internal complex type value (non-MP)
using MHDComplex = MHDComplex;
/// @brief Typedef for internal extended floating point type (non-MP)
using MHDLong = MHDLong;
#endif

///
/// START DEPRECATED use literal _mp instead
///
#ifdef QUICC_MULTPRECISION
/// Define a small macro to replace float constants to strings in the case of MP
/// computations
#define MHD_MP(c) MP::MHDFloat(#c)

/// Define a small macro to replace float constants to strings in the case of MP
/// computations
#define MHD_MP_LONG(c) MP::MHDFloat(#c)
#else
/// For normal computations the macro does nothing
#define MHD_MP(c) c
/// For normal computations the macro does nothing
#define MHD_MP_LONG(c) c##L
#endif // QUICC_MULTPRECISION
///
/// END DEPRECATED use literal _mp instead
///

} // namespace Internal
} // namespace QuICC

#endif // QUICC_TYPES_INTERNAL_BASICTYPES_HPP
