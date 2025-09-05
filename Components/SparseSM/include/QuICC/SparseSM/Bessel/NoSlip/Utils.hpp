/**
 * @file Utils.hpp
 * @brief Utils
 */

#ifndef QUICC_SPARSESM_BESSEL_NOSLIP_UTILS_HPP
#define QUICC_SPARSESM_BESSEL_NOSLIP_UTILS_HPP

// System includes
//

// Project includes
//
#include "Types/Internal/Typedefs.hpp"

namespace QuICC {

namespace SparseSM {

namespace Bessel {

namespace NoSlip {

/**
 * @brief Compute Bessel roots for no-slip boundary condition
 *
 * @param roots   Output vector of computed roots
 * @param l       Harmonic degree l
 * @param nRoots  Number of roots to compute
 */
void getRoots(std::vector<Internal::MHDFloat>& roots, const int l,
   const int nRoots);


} // namespace NoSlip
} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_NOSLIP_UTILS_HPP
