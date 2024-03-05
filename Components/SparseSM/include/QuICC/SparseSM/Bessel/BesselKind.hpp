/**
 * @file BesselKind.hpp
 * @brief Existing Bessel basis kinds
 */

#ifndef QUICC_SPARSESM_BESSEL_BESSELKIND_HPP
#define QUICC_SPARSESM_BESSEL_BESSELKIND_HPP

// Project includes
//

namespace QuICC {

namespace SparseSM {

namespace Bessel {

/// Different kinds of Bessel basis
enum class BesselKind
{
   VALUE,
   INSULATING,
};

} // namespace Bessel
} // namespace SparseSM
} // namespace QuICC

#endif // QUICC_SPARSESM_BESSEL_BESSELKIND_HPP
