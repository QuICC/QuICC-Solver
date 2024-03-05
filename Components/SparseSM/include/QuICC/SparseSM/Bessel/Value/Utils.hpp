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
 *
 * @param roots   Output vector of computed roots
 * @param l       Harmonic degree l
 * @param nRoots  Number of roots to compute
 */
void getRoots(std::vector<Internal::MHDFloat>& roots, const int l,
   const int nRoots);


} // namespace Value
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_VALUE_UTILS_HPP
