/**
 * @file Utils.hpp
 * @brief Utils
 */

#ifndef QUICC_SPARSESM_BESSEL_VALUE_UTILS_HPP
#define QUICC_SPARSESM_BESSEL_VALUE_UTILS_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace Value {

   /**
    * @brief Compute Bessel roots for zero value boundary condition
    */
   void getRoots(std::vector<Internal::MHDFloat>& roots, const int l, const int nRoots);


}
}
}
}

#endif // QUICC_SPARSESM_BESSEL_VALUE_UTILS_HPP
