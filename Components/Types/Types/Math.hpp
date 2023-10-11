/**
 * @file Math.hpp
 * @brief Definition of some useful math constants
 */

#ifndef QUICC_TYPES_MATH_HPP
#define QUICC_TYPES_MATH_HPP

// System includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Math {

/**
 * @brief The constant \f$\pi\f$
 */
constexpr MHDFloat PI = 3.14159265358979323846;

/**
 * @brief The constant \f$\pi\f$ for extended precision
 */
constexpr MHDLong PI_long = 3.14159265358979323846264338327950288419717l;

/**
 * @brief Pure imaginary value I
 */
constexpr MHDComplex cI{0.0, 1.0};

} // namespace Math
} // namespace QuICC

#endif // QUICC_TYPES_MATH_HPP
