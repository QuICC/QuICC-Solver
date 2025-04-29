/**
 * @file BasicTypes.hpp
 * @brief Some basic typedefs used in the whole project
 */

#ifndef QUICC_TYPES_BASICTYPES_HPP
#define QUICC_TYPES_BASICTYPES_HPP

// System includes
//
#include <complex>

// Project includes
//

namespace QuICC {

/// @brief Typedef for floating point type scalar
typedef double MHDFloat;
/// @brief Typedef for complex type scalar
typedef std::complex<MHDFloat> MHDComplex;
/// @brief Typedef for extended native floating point scalar
typedef long double MHDLong;

} // namespace QuICC

#endif // QUICC_BASICTYPES_HPP
