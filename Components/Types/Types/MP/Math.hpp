/**
 * @file Math.hpp
 * @brief Definition of multiple precision constants and math operations
 */

#ifndef QUICC_TYPES_MP_MATH_HPP
#define QUICC_TYPES_MP_MATH_HPP

// System includes
//

// Project includes
//
#include "Types/MP/BasicTypes.hpp"
#include "Types/MP/Literals.hpp"

namespace QuICC {
namespace MP {
/// @brief namespace for constants and math operations
namespace Math {

// bring in math ops
using namespace boost::multiprecision;
using namespace MP::Literals;

static_assert(QUICC_MULTPRECISION_DIGITS < 286, "PI literal is not accurate");

/// @brief The constant \f$\pi\f$ to 286 digits
const MHDFloat PI =
   3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856_mp;
/// @brief The constant \f$\pi\f$ to 286 digits
const MHDFloat PI_long = PI;
/// @brief Pure imaginary value I
const MHDComplex cI{0.0_mp, 1.0_mp};

} // namespace Math
} // namespace MP
} // namespace QuICC

#endif // QUICC_TYPES_MP_MATH_HPP
