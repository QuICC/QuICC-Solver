/**
 * @file BasicTypes.hpp
 * @brief Definitio
 *
 * /// @brief Complex typedef for MPn of multiple precision basic types
 */
/// @brief Complex typedef for MP

#ifndef QUICC_TYPES_MP_BASICTYPES_HPP
#define QUICC_TYPES_MP_BASICTYPES_HPP

// System includes
//
#include <complex>
#if defined QUICC_MPBACKEND_BOOST
#include <boost/multiprecision/cpp_dec_float.hpp>
#elif defined QUICC_MPBACKEND_GMP
#include <boost/multiprecision/gmp.hpp>
#elif defined QUICC_MPBACKEND_MPFR
#include <boost/multiprecision/mpfr.hpp>
#elif defined QUICC_MPBACKEND_QUAD
#include <boost/multiprecision/float128.hpp>
#endif // defined QUICC_MPBACKEND_BOOST
#include <boost/multiprecision/eigen.hpp>

// Project includes
//

namespace QuICC {
/// @brief namespace for multiple precision types and operations
namespace MP {

/// @brief Floating point typedef for MP
#if defined QUICC_MPBACKEND_BOOST
typedef boost::multiprecision::number<
   boost::multiprecision::cpp_dec_float<QUICC_MULTPRECISION_DIGITS>>
   MHDFloat;
#elif defined QUICC_MPBACKEND_GMP
typedef boost::multiprecision::number<
   boost::multiprecision::gmp_float<QUICC_MULTPRECISION_DIGITS>>
   MHDFloat;
#elif defined QUICC_MPBACKEND_MPFR
typedef boost::multiprecision::number<
   boost::multiprecision::mpfr_float_backend<QUICC_MULTPRECISION_DIGITS>>
   MHDFloat;
#elif defined QUICC_MPBACKEND_QUAD
typedef boost::multiprecision::float128 MHDFloat;
#endif

/// @brief Typedef for complex type value
typedef std::complex<MHDFloat> MHDComplex;

/// @brief Typedef for extended native floating point value
typedef MHDFloat MHDLong;

} // namespace MP
} // namespace QuICC

#endif // QUICC_TYPES_MP_BASICTYPES_HPP
