/**
 * @file Precision.cpp
 * @brief Implementation of precision related constants and routines
 */

// System includes
//

// External includes
//

// Class include
//
#include "Types/Precision.hpp"

// Project includes
//
#include "Types/Constants.hpp"

namespace QuICC {

   #ifdef QUICC_MULTPRECISION
      const internal::MHDFloat Precision::PI = precision::acos(MHD_MP(-1.0));

      const internal::MHDLong Precision::PI_long = precision::acos(MHD_MP_LONG(-1.0));

      const bool Precision::hasMP = true;
   #else
      const internal::MHDFloat Precision::PI = Math::PI;

      const internal::MHDLong Precision::PI_long = 3.14159265358979323846264338327950288419717l;

      const bool Precision::hasMP = false;
   #endif// QUICC_MULTPRECISION

}
