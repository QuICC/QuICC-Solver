/**
 * @file Math.hpp
 * @brief Definition of internal constants and math operations
 */

#ifndef QUICC_TYPES_INTERNAL_MATH_HPP
#define QUICC_TYPES_INTERNAL_MATH_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/BasicTypes.hpp"
#ifdef QUICC_MULTPRECISION
#include "Types/MP/Math.hpp"
#else
#include "Types/Math.hpp"
#endif

namespace QuICC {
namespace Internal {
/// @brief namespace for constants and math operations
namespace Math {

#ifdef QUICC_MULTPRECISION
// bring in math ops
using namespace MP::Math;
/// @brief The constant \f$\pi\f$
const Internal::MHDFloat PI = MP::Math::PI;
/// @brief The constant \f$\pi\f$ with extended precision
const Internal::MHDLong PI_long = MP::Math::PI_long;
/// @brief Pure imaginary value I
const Internal::MHDComplex cI = MP::Math::cI;
#else
// bring in math ops
using namespace std;
/// @brief The constant \f$\pi\f$
const Internal::MHDFloat PI = QuICC::Math::PI;
/// @brief The constant \f$\pi\f$ with extended precision
const Internal::MHDLong PI_long = QuICC::Math::PI_long;
/// @brief Pure imaginary value I
const Internal::MHDComplex cI = QuICC::Math::cI;
#endif

} // namespace Math
} // namespace Internal
} // namespace QuICC

#endif // QUICC_TYPES_INTERNAL_MATH_HPP
