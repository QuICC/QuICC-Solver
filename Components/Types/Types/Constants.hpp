/**
 * @file Constants.hpp
 * @brief Definition of some useful math constants
 */

#ifndef QUICC_MATH_CONSTANTS_HPP
#define QUICC_MATH_CONSTANTS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Types/Typedefs.hpp"

namespace QuICC {

namespace Math {

   /**
    * @brief The constant \f$\pi\f$
    */
   const MHDFloat PI = 3.14159265358979323846;

   /**
    * @brief Pure imaginary value I
    */
   const MHDComplex cI = MHDComplex(0.0, 1.0);

}
}

#endif // QUICC_MATH_CONSTANTS_HPP
