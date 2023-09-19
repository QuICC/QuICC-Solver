/**
 * @file BasicTypes.hpp
 * @brief Some basic typedefs used in the whole project 
 */

#ifndef QUICC_BASICTYPES_HPP
#define QUICC_BASICTYPES_HPP

// Configuration includes
//

// System includes
//
#include <complex>
#include <variant>

// External includes
//

// Project includes
//

namespace QuICC {

   /**
    * @name Basic scalar types typedefs
    */
   //@{
   /// Typedef for floating point type value
   typedef double MHDFloat;
   /// Typedef for complex type value
   typedef std::complex<MHDFloat> MHDComplex;
   /// Typedef for extended native floating point value
   typedef long double MHDLong;

#if defined __NVCC__
//#warning "Variant type is not supported until CMake 3.18"
#else
   /// Typedef for real/complex variant
   typedef std::variant<MHDFloat,MHDComplex> MHDVariant;
#endif
   //@}

}

#endif // QUICC_BASICTYPES_HPP
